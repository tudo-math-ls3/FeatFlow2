!##############################################################################
!# ****************************************************************************
!# <name> spacetimelinearsystem </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains routines to support a space-time linear system
!# for the underlying problem. As this is problem specific, there is a special
!# routine included here for setting up the and doing multiplication with
!# the global matrix.
!#
!# The routines in this module form the basic set of routines:
!#
!# 1.) cc_spaceTimeMatVec
!#     -> Matrix-Vector multiplication of a space-time vector with the
!#        global matrix.
!#
!# 2.) cc_generateInitCondRHS
!#     -> From the initial solution vector at time t=0, generate a
!#        RHS vector for the initial condition.
!#
!# Auxiliary routines:
!#
!# 1.) cc_setupMatrixWeights
!#     -> Initialise the matrix weights of a submatrix in the global space-time
!#        matrix.
!#
!# 2.) cc_disableSubmatrix
!#     -> Disables a submatrix by setting the corresponding weights to zero.
!#
!# 3.) cc_getFullMatrixDummy
!#     -> Get a dummy structure of a full matrix
!#
!# </purpose>
!#
!# The structure t_ccoptSpaceTimeDiscretisation realises in conjunction with
!# the above assembling routines a large space-time coupled system which solves
!# the nonstationary (Navier-)Stokes problem in all time steps simultaneously.
!# The main system matrix is assembled from the matrix blocks that are used in
!# each time step of a usual time-dependent simulation. It has the following
!# form (here for 2 timesteps with solutions y_0(initial), y_1 and y_2):
!#
!#  Explicit Euler.
!#  a=alpha. g=gamma. y_i=velocity in the i'th step. l_i=dual velocity in the i'th step.
!#
!#  Stokes:          A = A(y_i) = -nu*Laplace(.)
!#                   R = -M (from the RHS)
!#
!#  Navier-Stokes:   Standard iteration:
!#                   A = A(y_i) = -nu*Laplace(.) + y_i*grad(.)
!#                   N = N(y_i) = -nu*Laplace(.) - y_i*grad(.) + (.)*grad(y_i)
!#                   R = R(l_i) = -M
!#                   P = 1/a M
!#
!#                   Newton iteration:
!#                   A = A(y_i) = -nu*Laplace(.) + y_i*grad(.) + grad(y_i)*(.)
!#                   N = N(y_i) = -nu*Laplace(.) - y_i*grad(.) + grad(y_i)*(.)
!#                   R = R(l_i) = -M + l_i*grad(.) - grad(l_i)*(.)
!#                   P = P[a,b](l_i) = M~(l_i)
!#
!#        with M~(l_i)_kl = int (c(l_i) phi_l phi_k) dx
!#        and c(l_i) = 1 if a < -1/alpha l_i(x) < b and c(l_i) = 0 otherwise
!#        (derivative of the projection operator for constrains on the control)
!#
!#  ====================================================================================================
!#  [I     ]                       |                                  |
!#                                 |                                  |
!#            [I ]                 |                                  |
!#                                 |                                  |
!#                 [M/dt + N] [-B] |                 [-M/dt    ]      |
!#                                 |                                  |
!#                 [-B^t    ]      |                                  |
!#                                 |                                  |
!#  -------------------------------+----------------------------------+---------------------------------
!#  [-M/dt ]                       | [M/dt + A] [-B] [ P       ]      |
!#                                 |                                  |
!#                                 | [-B^t    ]                       |
!#                                 |                                  |
!#                                 | [ R      ]      [M/dt + N ] [-B] |                  [-M/dt   ]
!#                                 |   ^                              |
!#                                 | from the        [-B^t     ]      |
!#                                 | RHS                              |
!#  -------------------------------+----------------------------------+---------------------------------
!#                                 | [-M/dt   ]                       |  [M/dt + A] [-B] [ P      ]
!#                                 |                                  |
!#                                 |                                  |  [-B^t    ]
!#                                 |                                  |
!#                                 |                                  |  [ g R    ]      [M/dt + N] [-B]
!#                                 |                                  |    ^
!#                                 |                                  |  from the        [-B^t    ]
!#                                 |                                  |  RHS
!#  ====================================================================================================
!#
!# More general:
!#
!#  Crank Nicolson.
!#  a=alpha. g=gamma. y_i=(primal) velocity in the i'th step. T=theta=0.5
!#  Stokes:          A(y_i) = -nu*Laplace(.)
!#                   N(y_i) = not present
!#  Navier-Stokes:   A(y_i) = -nu*Laplace(.) + y_i*grad(.)
!#                   N(y_i) = -nu*Laplace(.) - y_i*grad(.) + (.)*grad(y_i)
!#
!# WARNING: WRONG. HAS TO BE REWRITTEN!
!#
!#  ===============================================================================================================================================================
!#  [I                ]                                   |                                                    |
!#                                                        |                                                    |
!#                      [I      ]                         |                                                    |
!#                                                        |                                                    |
!#                                [M/dt + T A(y_0)] [-B ] |                          [-M/dt+(1-T)N(y_0)]       |
!#                                                        |                                        ^           |
!#                                [-B^t           ]       |                                        |           |
!#                                                        |                                    !correct!       |
!#  ------------------------------------------------------+----------------------------------------------------+---------------------------------------------------
!#  [-M/dt+(1-T)A(y_0)]                                   | [M/dt + T A(y_1) ] [-B ] [ P               ]       |
!#                                                        |                                                    |
!#                                                        | [-B^t            ]                                 |
!#                                                        |                                                    |
!#                                                        | [-T M            ]       [ M/dt + T N(y_1) ] [-B ] |                          [-M/dt+(1-T)N(y_1)]
!#                                                        |   ^                                                |                                        ^
!#                                                        | from the                 [ -B^t            ]       |                                        |
!#                                                        | RHS                                                |                                    !correct!
!#  ------------------------------------------------------+----------------------------------------------------+---------------------------------------------------
!#                                                        | [-M/dt+(1-T)A(y_1)]                                |  [M/dt + T A(y_2)] [-B ] [ P               ]
!#                                                        |                                                    |
!#                                                        |                                                    |  [-B^t           ]
!#                                                        |                                                    |
!#                                                        |                                                    |  [-g T M         ]       [ M/dt + T N(y_2) ] [-B ]
!#                                                        |                                                    |    ^
!#                                                        |                                                    |  from the                [ -B^t            ]
!#                                                        |                                                    |  RHS
!#  ===============================================================================================================================================================
!#
!##############################################################################

module spacetimelinearsystem

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
  use linearsolverautoinitialise
  use matrixrestriction
  use paramlist
  use timestepping
  
  use collection
  use convection
    
  use basicstructures
  use user_callback

  !use spacepreconditioner
  !use spacepreconditionerinit
  use timeanalysis
  use spatialbc
  use spacediscretisation
  use spacematvecassembly
  
  use timediscretisation
  use spacetimediscretisation
  use spacetimevectors
  use timeboundaryconditions
  
  use matrixio
    
  implicit none

!<types>

!<typeblock>

  ! Defines the basic shape of the supersystem which realises the coupling
  ! between all timesteps.
  type t_ccoptSpaceTimeMatrix
  
    ! Type of the matrix.
    ! =0: Standard system matrix.
    ! =1: Newton matrix
    integer :: cmatrixType = 0
    
    ! Activate nonlinear projection operators in the matrix.
    ! =0: Matrix does not contain projection operators for u.
    ! =1: Matrix contains projection operators for u.
    integer :: ccontrolConstraints = 0
  
    ! Pointer to the space time discretisation that corresponds to this matrix.
    type(t_ccoptSpaceTimeDiscretisation), pointer :: p_rspaceTimeDiscr => null()
  
    ! Pointer to a space-time solution vector that defines the point
    ! where the nonlinearity is evaluated when applying or calculating
    ! matrices.
    ! If this points to NULL(), there is no global solution specified,
    ! so the matrix assembly/application routines use the vector that
    ! is specified in their input parameters.
    type(t_spacetimeVector), pointer :: p_rsolution => null()
    
  end type

!</typeblock>

!</types>

!<constants>

!<constantblock description="Filter-ID's of the filters to be applied in a matrix vector multiplication">

  ! No filter
  integer(I32), parameter :: SPTID_FILTER_NONE  = 0

  ! Filter the output vector for boundary conditions in the defect.
  integer(I32), parameter :: SPTID_FILTER_BCDEF = 2**0
  
  ! Filter the output vector for initial conditions in the defect
  integer(I32), parameter :: SPTID_FILTER_ICDEF = 2**1

  ! Filter the output vector for terminal conditions in the defect
  integer(I32), parameter :: SPTID_FILTER_TCDEF = 2**2
  
  ! Apply all filters that are typically applied to a defect vector
  integer(I32), parameter :: SPTID_FILTER_DEFECT = SPTID_FILTER_BCDEF + &
                                                   SPTID_FILTER_ICDEF + &
                                                   SPTID_FILTER_TCDEF
                                                   
!</constantblock>

!</constants>

contains

  
  ! ***************************************************************************
  
!<subroutine>

  subroutine cc_setupMatrixWeights (rproblem,rspaceTimeMatrix,dtheta,&
      isubstep,irelpos,rnonlinearSpatialMatrix)

!<description>
  ! This routine sets up the matrix weights in the rnonlinearSpatialMatrix structure
  ! according to the position of the corresponding submatrix.
  ! isubstep defines the timestep which is to be tackled -- which corresponds
  ! to the row in the supermatrix.
  ! irelpos specifies the column in the supermatrix relative to the diagonal --
  ! thus =0 means the diagonal, =1 the submatrix above the diagonal and =-1
  ! the submatrix below the diagonal.
  !
  ! The matrix weights in rnonlinearSpatialMatrix are initialised based on this
  ! information such that the 'assembleMatrix' and 'assembleDefect' routines of
  ! the core equation will build the system matrix / defect at that position
  ! in the supermatrix.
  !
  ! The routine does not initialise the pointers to the basic matrices/discretisation
  ! structures in the rnonlinearSpatialMatrix structure. This has to be done by
  ! the caller!
!</description>

!<input>
  ! Problem structure
  type(t_problem), intent(INOUT) :: rproblem
  
  ! A t_ccoptSpaceTimeMatrix structure defining the discretisation of the
  ! coupled space-time matrix.
  type(t_ccoptSpaceTimeMatrix), intent(IN), target :: rspaceTimeMatrix

  ! Theta scheme identifier.
  ! = 0.5: Crank-Nicolson.
  ! = 1.0: Explicit Euler
  real(DP), intent(IN) :: dtheta
  
  ! Substep in the time-dependent simulation = row in the supermatrix.
  ! Range 0..nsubsteps
  integer, intent(IN) :: isubstep
  
  ! Specifies the column in the supermatrix relative to the diagonal.
  ! =0: set matrix weights for the diagonal.
  integer, intent(IN) :: irelpos
!</input>

!<inputoutput>
  ! A t_nonlinearSpatialMatrix structure that defines the shape of the core
  ! equation. The weights that specify the submatrices of a small 6x6
  ! block matrix system are initialised depending on the position
  ! specified by isubstep and nsubsteps.
  !
  ! The structure must have been initialised with cc_initNonlinMatrix!
  type(t_nonlinearSpatialMatrix), intent(inout) :: rnonlinearSpatialMatrix
!</inputoutput>

!</subroutine>

!    ! If the following constant is set from 1.0 to 0.0, the primal system is
!    ! decoupled from the dual system!
!    real(DP), parameter :: dprimalDualCoupling = 1.0_DP
!
!    ! If the following constant is set from 1.0 to 0.0, the dual system is
!    ! decoupled from the primal system!
!    real(DP), parameter :: ddualPrimalCoupling = 1.0_DP
!
!    ! If the following parameter is set from 1.0 to 0.0, the terminal
!    ! condition between the primal and dual equation is decoupled, i.e.
!    ! the dual equation gets independent from the primal one.
!    real(DP), parameter :: dterminalCondDecoupled = 1.0_DP
!
!    ! If the following parameter is set from 1.0 to 0.0, the time coupling
!    ! is disabled, resulting in a stationary simulation in every timestep.
!    real(DP), parameter :: dtimeCoupling = 1.0_DP

    real(DP) :: dprimalDualCoupling,ddualPrimalCoupling
    real(DP) :: dterminalCondDecoupled,dtimeCoupling

    ! This constant defines the type of equation. There are two equivalent
    ! formulations of the dual equation which only differs in the sign
    ! of the dual velocity.
    ! A constant of "1" here means to use the formulation with "y-z" in
    ! the RHS of the dual equation, while a constant of "-1" means to use the
    ! formulation with "-(y-z)" there.
    ! Note: If this is changed, a "-" sign must be implemented / removed from
    ! the RHS, too!
    real(DP) :: dequationType
    
    ! Pointer to the space time discretisation structure.
    type(t_ccoptSpaceTimeDiscretisation), pointer :: p_rspaceTimeDiscr
    real(DP) :: dnewton
    real(DP) :: dtstep
    logical :: bconvectionExplicit
    
    call parlst_getvalue_double (rproblem%rparamList,'DEBUG',&
        'dprimalDualCoupling',dprimalDualCoupling,1.0_DP)

    call parlst_getvalue_double (rproblem%rparamList,'DEBUG',&
        'ddualPrimalCoupling',ddualPrimalCoupling,1.0_DP)

    call parlst_getvalue_double (rproblem%rparamList,'DEBUG',&
        'dterminalCondDecoupled',dterminalCondDecoupled,1.0_DP)

    call parlst_getvalue_double (rproblem%rparamList,'DEBUG',&
        'dtimeCoupling',dtimeCoupling,1.0_DP)

    p_rspaceTimeDiscr => rspaceTimeMatrix%p_rspaceTimeDiscr
    
    dequationType = 1.0_DP
    if (rproblem%roptcontrol%ispaceTimeFormulation .ne. 0) &
      dequationType = -1.0_DP
      
    ! Treat the convection explicitely?
    bconvectionExplicit = rproblem%roptcontrol%iconvectionExplicit .ne. 0

    ! What's the matrix type we should set up? If we have to set up a
    ! Newton matrix, we put dnewton to 1.0 so the Newton part in the
    ! primal velocity is assembled.
    dnewton = 0.0_DP
    if (rspaceTimeMatrix%cmatrixType .eq. 1) then
      if (rproblem%rphysicsPrimal%iequation .eq. 0) then
        ! Newton is only to be assembled in Navier-Stokes!
        dnewton = 1.0_DP
      end if
    end if
    
    dtstep = p_rspaceTimeDiscr%rtimeDiscr%dtstep

    ! Clear the coefficients
    rnonlinearSpatialMatrix%Diota(:,:) = 0.0_DP
    rnonlinearSpatialMatrix%Dalpha(:,:) = 0.0_DP
    rnonlinearSpatialMatrix%Dtheta(:,:) = 0.0_DP
    rnonlinearSpatialMatrix%Dgamma(:,:) = 0.0_DP
    rnonlinearSpatialMatrix%Dnewton(:,:) = 0.0_DP
    rnonlinearSpatialMatrix%DgammaT(:,:) = 0.0_DP
    rnonlinearSpatialMatrix%Dnewton2(:,:) = 0.0_DP
    rnonlinearSpatialMatrix%DgammaT2(:,:) = 0.0_DP
    rnonlinearSpatialMatrix%DnewtonT(:,:) = 0.0_DP
    rnonlinearSpatialMatrix%Deta(:,:) = 0.0_DP
    rnonlinearSpatialMatrix%Dtau(:,:) = 0.0_DP
    rnonlinearSpatialMatrix%Dkappa(:,:) = 0.0_DP

    ! The first and last substep is a little bit special concerning
    ! the matrix!
    if (isubstep .eq. 0) then
      
      ! We are in the first substep
      !
      ! Specify which vectors should be used for evaluating the nonlinearity.
      ! Aii-1 = undefined.
      ! Aii = Aii(x_i,lambda_i).
      ! Aii+1 = Aii+1(x_i,lambda_i+1)
    
      if (irelpos .eq. 0) then
      
        ! The diagonal matrix.
        if (.not. bconvectionExplicit) then
          rnonlinearSpatialMatrix%iprimalSol = 2
          rnonlinearSpatialMatrix%idualSol = 2
          rnonlinearSpatialMatrix%idualSol2 = 3
        else
          rnonlinearSpatialMatrix%iprimalSol = 2
          rnonlinearSpatialMatrix%idualSol = 3
          rnonlinearSpatialMatrix%idualSol2 = 2
        end if

        rnonlinearSpatialMatrix%Dalpha(1,1) = dtimeCoupling * 1.0_DP/dtstep
        rnonlinearSpatialMatrix%Dalpha(2,2) = dtimeCoupling * 1.0_DP/dtstep
        
        rnonlinearSpatialMatrix%Dtheta(1,1) = 1.0_DP
        rnonlinearSpatialMatrix%Dtheta(2,2) = 1.0_DP
        
        if (.not. bconvectionExplicit) then

          rnonlinearSpatialMatrix%Dgamma(1,1) = &
              real(1-rproblem%rphysicsPrimal%iequation,DP)
          rnonlinearSpatialMatrix%Dgamma(2,2) = &
              - real(1-rproblem%rphysicsPrimal%iequation,DP)
          
          rnonlinearSpatialMatrix%Dnewton(1,1) = dnewton
          rnonlinearSpatialMatrix%DnewtonT(2,2) = &
                real(1-rproblem%rphysicsPrimal%iequation,DP)
                
        end if

        rnonlinearSpatialMatrix%Deta(1,1) = 1.0_DP
        rnonlinearSpatialMatrix%Deta(2,2) = 1.0_DP
        
        rnonlinearSpatialMatrix%Dtau(1,1) = 1.0_DP
        rnonlinearSpatialMatrix%Dtau(2,2) = 1.0_DP
        
        ! Only difference to the usual diagonal: No coupling mass matrix
        ! from dual to the primal velocity.
        rnonlinearSpatialMatrix%Dalpha(2,1) = ddualPrimalCoupling * &
            (-dequationType)
        
        if (.not. bconvectionExplicit) then

          if (dnewton .ne. 0.0_DP) then
            rnonlinearSpatialMatrix%DgammaT(2,1) = ddualPrimalCoupling * &
                ( dequationType)
            !rnonlinearSpatialMatrix%Dgamma(2,1) = ddualPrimalCoupling * &
            !    ( dequationType)
            rnonlinearSpatialMatrix%Dnewton(2,1) = ddualPrimalCoupling * &
                (-dequationType)

            ! For Crank-Nicolson there appears a 2nd reactive term
            ! stemming from the next timestep.

            rnonlinearSpatialMatrix%DgammaT2(2,1) = ddualPrimalCoupling * &
                (1.0_DP-dtheta) * ( dequationType)
            !rnonlinearSpatialMatrix%Dgamma2(2,1) = ddualPrimalCoupling * &
            !    (1.0_DP-dtheta) * ( dequationType)
            rnonlinearSpatialMatrix%Dnewton2(2,1) = ddualPrimalCoupling * &
                (1.0_DP-dtheta) * (-dequationType)
          end if
          
        else
        
          if (dnewton .ne. 0.0_DP) then
            rnonlinearSpatialMatrix%DgammaT(2,1) = ddualPrimalCoupling * &
                ( dequationType)
            !rnonlinearSpatialMatrix%Dgamma(2,1) = ddualPrimalCoupling * &
            !    ( dequationType)
            rnonlinearSpatialMatrix%Dnewton(2,1) = ddualPrimalCoupling * &
                (-dequationType)
          end if
        
        end if
                  
      else if (irelpos .eq. 1) then
      
        ! Offdiagonal matrix on the right of the diagonal.
        !
        ! Note that at this point, the nonlinearity must be evaluated
        ! at xi due to the discretisation scheme!!!
        !
        
        !wirklich? unten nochmal!

        rnonlinearSpatialMatrix%iprimalSol = 2
        rnonlinearSpatialMatrix%idualSol = 2
        rnonlinearSpatialMatrix%idualSol2 = 3
        
        ! WARNING: For a very strange situation, taking xi here (which is said
        ! to be the correct evaluation point from the theory) does not lead
        ! to quadratic convergence in the Newton. Taking xi+1 does!?!
        ! So we take xi+1, although the theory tells us to take xi!
        ! ...
        ! No, that does not to be right. Commented out since the above works
        ! as well and should be correct due to the theory.
        ! rnonlinearSpatialMatrix%idualSol = 3
        
        ! Switch off any stabilisation
        rnonlinearSpatialMatrix%dupsam1 = 0.0_DP
        rnonlinearSpatialMatrix%dupsam2 = 0.0_DP

        ! Create the matrix
        !   -M + dt*dtheta*[-nu\Laplace u + u \grad u]

        rnonlinearSpatialMatrix%Dalpha(2,2) = dtimeCoupling * (-1.0_DP)/dtstep
        
        rnonlinearSpatialMatrix%Dtheta(2,2) = (1.0_DP-dtheta)
        
        if (.not. bconvectionExplicit) then
        
          rnonlinearSpatialMatrix%Dgamma(2,2) = &
              - (1.0_DP-dtheta) * real(1-rproblem%rphysicsPrimal%iequation,DP)
          
          rnonlinearSpatialMatrix%DnewtonT(2,2) = &
                (1.0_DP-dtheta) * real(1-rproblem%rphysicsPrimal%iequation,DP)
        else

          rnonlinearSpatialMatrix%Dgamma(2,2) = &
              - dtheta * real(1-rproblem%rphysicsPrimal%iequation,DP)
          
          rnonlinearSpatialMatrix%DnewtonT(2,2) = &
                dtheta * real(1-rproblem%rphysicsPrimal%iequation,DP)

        end if

        rnonlinearSpatialMatrix%Dalpha(2,1) = ddualPrimalCoupling * &
            (-dequationType) * (1.0_DP-dtheta)
            
      end if
    
    else if (isubstep .lt. p_rspaceTimeDiscr%NEQtime-1) then
      
      ! We are sonewhere in the middle of the matrix. There is a substep
      ! isubstep+1 and a substep isubstep-1!
      
      ! -----
      
      ! Specify which vectors should be used for evaluating the nonlinearity.
      ! Aii-1 = Aii+1(x_i-1,lambda_i-1).
      ! Aii = Aii(x_i,lambda_i).
      ! Aii+1 = Aii+1(x_i,lambda_i+1)

      if (irelpos .eq. -1) then
      
        ! Matrix on the left of the diagonal.

        rnonlinearSpatialMatrix%iprimalSol = 1
        rnonlinearSpatialMatrix%idualSol = 1

        ! Switch off any stabilisation
        rnonlinearSpatialMatrix%dupsam1 = 0.0_DP
        rnonlinearSpatialMatrix%dupsam2 = 0.0_DP

        ! Create the matrix
        !   -M + dt*dtheta*[-nu\Laplace u + u \grad u]

        rnonlinearSpatialMatrix%Dalpha(1,1) = dtimeCoupling * (-1.0_DP)/dtstep
        
        rnonlinearSpatialMatrix%Dtheta(1,1) = (1.0_DP-dtheta)
        
        if (.not. bconvectionExplicit) then

          rnonlinearSpatialMatrix%Dgamma(1,1) = &
              (1.0_DP-dtheta) * real(1-rproblem%rphysicsPrimal%iequation,DP)
          
          rnonlinearSpatialMatrix%Dnewton(1,1) = &
              (1.0_DP-dtheta) * dnewton

        else
        
          rnonlinearSpatialMatrix%Dgamma(1,1) = &
              dtheta * real(1-rproblem%rphysicsPrimal%iequation,DP)
          
          rnonlinearSpatialMatrix%Dnewton(1,1) = &
              dtheta * dnewton
        
        end if

        rnonlinearSpatialMatrix%Dalpha(1,2) = dprimalDualCoupling * &
            dequationType * (1.0_DP-dtheta) / p_rspaceTimeDiscr%dalphaC

      else if (irelpos .eq. 0) then

        ! The diagonal matrix.
        if (.not. bconvectionExplicit) then
          rnonlinearSpatialMatrix%iprimalSol = 2
          rnonlinearSpatialMatrix%idualSol = 2
          rnonlinearSpatialMatrix%idualSol2 = 3
        else
          rnonlinearSpatialMatrix%iprimalSol = 2
          rnonlinearSpatialMatrix%idualSol = 3
          rnonlinearSpatialMatrix%idualSol2 = 2
        end if

        rnonlinearSpatialMatrix%Dalpha(1,1) = dtimeCoupling * 1.0_DP/dtstep
        rnonlinearSpatialMatrix%Dalpha(2,2) = dtimeCoupling * 1.0_DP/dtstep
        
        rnonlinearSpatialMatrix%Dtheta(1,1) = dtheta
        rnonlinearSpatialMatrix%Dtheta(2,2) = dtheta
        
        if (.not. bconvectionExplicit) then

          rnonlinearSpatialMatrix%Dgamma(1,1) = &
              dtheta * real(1-rproblem%rphysicsPrimal%iequation,DP)
          rnonlinearSpatialMatrix%Dgamma(2,2) = &
              - dtheta * real(1-rproblem%rphysicsPrimal%iequation,DP)
          
          rnonlinearSpatialMatrix%Dnewton(1,1) = dtheta * dnewton
          rnonlinearSpatialMatrix%DnewtonT(2,2) = &
                dtheta * real(1-rproblem%rphysicsPrimal%iequation,DP)
                
        end if

        rnonlinearSpatialMatrix%Deta(1,1) = 1.0_DP
        rnonlinearSpatialMatrix%Deta(2,2) = 1.0_DP
        
        rnonlinearSpatialMatrix%Dtau(1,1) = 1.0_DP
        rnonlinearSpatialMatrix%Dtau(2,2) = 1.0_DP
        
        rnonlinearSpatialMatrix%Dalpha(1,2) = dprimalDualCoupling * &
            dequationType * dtheta * 1.0_DP / p_rspaceTimeDiscr%dalphaC
        rnonlinearSpatialMatrix%Dalpha(2,1) = ddualPrimalCoupling * &
            (-dequationType) * dtheta
            
        if (.not. bconvectionExplicit) then

          if (dnewton .ne. 0.0_DP) then
            rnonlinearSpatialMatrix%DgammaT(2,1) = ddualPrimalCoupling * &
                ( dequationType) * dtheta
            !rnonlinearSpatialMatrix%Dgamma(2,1) = ddualPrimalCoupling * &
            !    ( dequationType) * dtheta
            rnonlinearSpatialMatrix%Dnewton(2,1) = ddualPrimalCoupling * &
                (-dequationType) * dtheta

            ! For Crank-Nicolson there appears a 2nd reactive term
            ! stemming from the next timestep.

            rnonlinearSpatialMatrix%DgammaT2(2,1) = ddualPrimalCoupling * &
                (1.0_DP-dtheta) * ( dequationType)
            !rnonlinearSpatialMatrix%Dgamma2(2,1) = ddualPrimalCoupling * &
            !    (1.0_DP-dtheta) * ( dequationType)
            rnonlinearSpatialMatrix%Dnewton2(2,1) = ddualPrimalCoupling * &
                (1.0_DP-dtheta) * (-dequationType)
          end if
          
        else
        
          if (dnewton .ne. 0.0_DP) then
            rnonlinearSpatialMatrix%DgammaT(2,1) = ddualPrimalCoupling * &
                ( dequationType) * dtheta
            !rnonlinearSpatialMatrix%Dgamma(2,1) = ddualPrimalCoupling * &
            !    ( dequationType) * dtheta
            rnonlinearSpatialMatrix%Dnewton(2,1) = ddualPrimalCoupling * &
                (-dequationType) * dtheta
          end if
        
        end if

      else if (irelpos .eq. 1) then
            
        ! Matrix on the right of the diagonal.
        !
        ! Note that at this point, the nonlinearity must be evaluated
        ! at xi due to the discretisation scheme!!!
        rnonlinearSpatialMatrix%iprimalSol = 2
        rnonlinearSpatialMatrix%idualSol = 2
        rnonlinearSpatialMatrix%idualSol2 = 3
        
        ! WARNING: For a very strange situation, taking xi here (which is said
        ! to be the correct evaluation point from the theory) does not lead
        ! to quadratic convergence in the Newton. Taking xi+1 does!?!
        ! So we take xi+1, although the theory tells us to take xi!
        ! ...
        ! No, that does not to be right. Commented out since the above works
        ! as well and should be correct due to the theory.
        ! rnonlinearSpatialMatrix%idualSol = 3

        ! Switch off any stabilisation
        rnonlinearSpatialMatrix%dupsam1 = 0.0_DP
        rnonlinearSpatialMatrix%dupsam2 = 0.0_DP

        ! Create the matrix
        !   -M + dt*dtheta*[-nu\Laplace u + u \grad u]
        rnonlinearSpatialMatrix%Dalpha(2,2) = dtimeCoupling * (-1.0_DP)/dtstep
        
        rnonlinearSpatialMatrix%Dtheta(2,2) = (1.0_DP-dtheta)
        
        if (.not. bconvectionExplicit) then

          rnonlinearSpatialMatrix%Dgamma(2,2) = &
              - (1.0_DP-dtheta) * real(1-rproblem%rphysicsPrimal%iequation,DP)
          
          rnonlinearSpatialMatrix%DnewtonT(2,2) = &
              (1.0_DP-dtheta) * real(1-rproblem%rphysicsPrimal%iequation,DP)
              
        else
        
          rnonlinearSpatialMatrix%Dgamma(2,2) = &
              - dtheta * real(1-rproblem%rphysicsPrimal%iequation,DP)
          
          rnonlinearSpatialMatrix%DnewtonT(2,2) = &
                dtheta * real(1-rproblem%rphysicsPrimal%iequation,DP)
        
        end if

        rnonlinearSpatialMatrix%Dalpha(2,1) = ddualPrimalCoupling * &
            (-dequationType) * (1.0_DP-dtheta)
            
        if (bconvectionExplicit) then

          ! DON'T KNOW IF THIS IS CORRECT!!!
        
          if (dnewton .ne. 0.0_DP) then
            rnonlinearSpatialMatrix%DgammaT(2,1) = ddualPrimalCoupling * &
                (1.0_DP-dtheta) * ( dequationType)
            !rnonlinearSpatialMatrix%Dgamma(2,1) = ddualPrimalCoupling * &
            !    (1.0_DP-dtheta) * ( dequationType)
            rnonlinearSpatialMatrix%Dnewton(2,1) = ddualPrimalCoupling * &
                (1.0_DP-dtheta) * (-dequationType)
          end if
          
        end if
            
      end if
    
    else
    
      ! We are in the last substep
      
      ! Specify which vectors should be used for evaluating the nonlinearity.
      ! Aii-1 = Aii+1(x_i-1,lambda_i-1).
      ! Aii = Aii(x_i,lambda_i).
      ! Aii+1 = Aii+1(x_i,lambda_i+1)

      if (irelpos .eq. -1) then
      
        ! Matrix on the left of the diagonal.
        
        rnonlinearSpatialMatrix%iprimalSol = 1
        rnonlinearSpatialMatrix%idualSol = 1

        ! Switch off any stabilisation
        rnonlinearSpatialMatrix%dupsam1 = 0.0_DP
        rnonlinearSpatialMatrix%dupsam2 = 0.0_DP

        ! Create the matrix
        !   -M + dt*dtheta*[-nu\Laplace u + u \grad u]

        rnonlinearSpatialMatrix%Dalpha(1,1) = dtimeCoupling * (-1.0_DP)/dtstep
        
        rnonlinearSpatialMatrix%Dtheta(1,1) = (1.0_DP-dtheta)
        
        if (.not. bconvectionExplicit) then

          rnonlinearSpatialMatrix%Dgamma(1,1) = &
              (1.0_DP-dtheta) * real(1-rproblem%rphysicsPrimal%iequation,DP)
          
          rnonlinearSpatialMatrix%Dnewton(1,1) = (1.0_DP-dtheta) * dnewton
          
        else
        
          rnonlinearSpatialMatrix%Dgamma(1,1) = &
              dtheta * real(1-rproblem%rphysicsPrimal%iequation,DP)
          
          rnonlinearSpatialMatrix%Dnewton(1,1) = dtheta * dnewton
          
        end if

        rnonlinearSpatialMatrix%Dalpha(1,2) = dprimalDualCoupling * &
            dequationType * (1.0_DP-dtheta) / p_rspaceTimeDiscr%dalphaC

      else if (irelpos .eq. 0) then

        ! The diagonal matrix.
        rnonlinearSpatialMatrix%iprimalSol = 2
        rnonlinearSpatialMatrix%idualSol = 2
        rnonlinearSpatialMatrix%idualSol2 = 2

        ! Current formulation:
        ! -(gamma+1/dt)*M*y + (M+dt*nu*L)*lambda = -(gamma+1/dt)*z
        
        rnonlinearSpatialMatrix%Dalpha(1,1) = dtimeCoupling * 1.0_DP/dtstep
        rnonlinearSpatialMatrix%Dalpha(2,2) = 1.0_DP/dtstep
        
        ! No 'time coupling' here; because of the terminal condition,
        ! the mass matrix resembles not the time dependence!
        
        rnonlinearSpatialMatrix%Dtheta(1,1) = dtheta
        rnonlinearSpatialMatrix%Dtheta(2,2) = dtheta
        
        if (.not. bconvectionExplicit) then

          rnonlinearSpatialMatrix%Dgamma(1,1) = &
              dtheta * real(1-rproblem%rphysicsPrimal%iequation,DP)
          rnonlinearSpatialMatrix%Dgamma(2,2) = &
              - dtheta * real(1-rproblem%rphysicsPrimal%iequation,DP)
          
          rnonlinearSpatialMatrix%Dnewton(1,1) = dtheta * dnewton
          rnonlinearSpatialMatrix%DnewtonT(2,2) = &
                dtheta * real(1-rproblem%rphysicsPrimal%iequation,DP)
                
        end if

        rnonlinearSpatialMatrix%Deta(1,1) = 1.0_DP
        rnonlinearSpatialMatrix%Deta(2,2) = 1.0_DP
        
        rnonlinearSpatialMatrix%Dtau(1,1) = 1.0_DP
        rnonlinearSpatialMatrix%Dtau(2,2) = 1.0_DP
        
        rnonlinearSpatialMatrix%Dalpha(1,2) = dprimalDualCoupling * &
            dequationType * dtheta * 1.0_DP / p_rspaceTimeDiscr%dalphaC
        rnonlinearSpatialMatrix%Dalpha(2,1) = ddualPrimalCoupling * &
            (-dequationType) * (dtheta + p_rspaceTimeDiscr%dgammaC / dtstep)
            
        if (.not. bconvectionExplicit) then
        
          ! Weight the mass matrix by GAMMA instead of delta(T).
          ! That's the only difference to the implementation above!
          if (dnewton .ne. 0.0_DP) then
            rnonlinearSpatialMatrix%DgammaT(2,1) = ddualPrimalCoupling * &
                ( dequationType) * dtheta
            !rnonlinearSpatialMatrix%Dgamma(2,1) = ddualPrimalCoupling * &
            !    ( dequationType) * dtheta
            rnonlinearSpatialMatrix%Dnewton(2,1) = ddualPrimalCoupling * &
                (-dequationType) * dtheta
          end if

        end if

      end if
        
    end if
    
    ! The dumin/dumax parameters are the same for all equations.
    rnonlinearSpatialMatrix%dumin1 = rproblem%roptcontrol%dumin1
    rnonlinearSpatialMatrix%dumax1 = rproblem%roptcontrol%dumax1
    rnonlinearSpatialMatrix%dumin2 = rproblem%roptcontrol%dumin2
    rnonlinearSpatialMatrix%dumax2 = rproblem%roptcontrol%dumax2
    
    ! General parameters.
    rnonlinearSpatialMatrix%dalphaC = rproblem%roptcontrol%dalphaC
    rnonlinearSpatialMatrix%ccontrolConstraints = rspaceTimeMatrix%ccontrolConstraints
    rnonlinearSpatialMatrix%cmatrixType = rspaceTimeMatrix%cmatrixType

  end subroutine
  
  ! ***************************************************************************
  
!<subroutine>

  subroutine cc_disableSubmatrix (rnonlinearSpatialMatrix,irow,icolumn)

!<description>
  ! Disables a subbklock in the nonlinear matrix rnonlinearSpatialMatrix.
  ! All weights of the correspopnding subblock are set to 0.
!</description>

!<input>
  ! The row/column of the submatrix to be disabled.
  integer :: irow,icolumn
!</input>

!<inputoutput>
  ! A t_nonlinearSpatialMatrix structure that defines the shape of the core
  ! equation. The weights that specify the submatrices of a small 6x6
  ! block matrix system are initialised depending on the position
  ! specified by isubstep and nsubsteps.
  !
  ! The structure must have been initialised with cc_initNonlinMatrix!
  type(t_nonlinearSpatialMatrix), intent(inout) :: rnonlinearSpatialMatrix
!</inputoutput>

!</subroutine>

    ! Clear the coefficients
    rnonlinearSpatialMatrix%Diota(irow,icolumn) = 0.0_DP
    rnonlinearSpatialMatrix%Dalpha(irow,icolumn) = 0.0_DP
    rnonlinearSpatialMatrix%Dtheta(irow,icolumn) = 0.0_DP
    rnonlinearSpatialMatrix%Dgamma(irow,icolumn) = 0.0_DP
    rnonlinearSpatialMatrix%Dnewton(irow,icolumn) = 0.0_DP
    rnonlinearSpatialMatrix%DgammaT(irow,icolumn) = 0.0_DP
    rnonlinearSpatialMatrix%Dnewton2(irow,icolumn) = 0.0_DP
    rnonlinearSpatialMatrix%DgammaT2(irow,icolumn) = 0.0_DP
    rnonlinearSpatialMatrix%DnewtonT(irow,icolumn) = 0.0_DP
    rnonlinearSpatialMatrix%Deta(irow,icolumn) = 0.0_DP
    rnonlinearSpatialMatrix%Dtau(irow,icolumn) = 0.0_DP
    rnonlinearSpatialMatrix%Dkappa(irow,icolumn) = 0.0_DP

  end subroutine

  ! ***************************************************************************

!<subroutine>

  subroutine cc_getFullMatrixDummy(rphysicsPrimal,rnonlinearSpatialMatrix)
  
!<description>
  ! Configures the parameters in rnonlinearSpatialMatrix with dummy
  ! coefficients so that the configuration resembles the matrix structure
  ! in the current situation. This is used during the allocation of memory
  ! to figure out the matrix structure and allocate necessary submatrices.
!</description>

!<input>
  ! Current physics structure of the primal equation
  type(t_physicsPrimal), intent(in) :: rphysicsPrimal
!</input>

!<inputoutput>
  ! A t_nonlinearSpatialMatrix structure that defines the shape of the core
  ! equation. The weights that specify the submatrices of a small 6x6
  ! block matrix system are initialised depending on the position
  ! specified by isubstep and nsubsteps. The coefficients in this
  ! structure are either set to 1.0 or 0.0 depending on whether a
  ! submatrix is active or not.
  !
  ! The structure must have been initialised with cc_initNonlinMatrix!
  type(t_nonlinearSpatialMatrix), intent(inout) :: rnonlinearSpatialMatrix
!</inputoutput>

!</subroutine>

    rnonlinearSpatialMatrix%Dtheta(1,1) = 1.0_DP   ! A velocity block
    rnonlinearSpatialMatrix%Deta(1,1) = 1.0_DP     ! A gradient block
    rnonlinearSpatialMatrix%Dtau(1,1) = 1.0_DP     ! A divergence block
    rnonlinearSpatialMatrix%Dkappa(1,1) = 1.0_DP   ! Pressure block
    ! A Newton block, if we have Navier-Stokes. For the case
    ! that we use Newton
    rnonlinearSpatialMatrix%Dnewton(1,1) = real(1-rphysicsPrimal%iequation,DP)
    
    rnonlinearSpatialMatrix%Dtheta(2,2) = 1.0_DP   ! A velocity block
    rnonlinearSpatialMatrix%Deta(2,2) = 1.0_DP     ! A gradient block
    rnonlinearSpatialMatrix%Dtau(2,2) = 1.0_DP     ! A divergence block
    ! A Newton block, if we have Navier-Stokes
    rnonlinearSpatialMatrix%Dnewton(2,2) = real(1-rphysicsPrimal%iequation,DP)
    rnonlinearSpatialMatrix%Dkappa(2,2) = 1.0_DP   ! Pressure block
    
    rnonlinearSpatialMatrix%Dalpha(1,1) = 1.0_DP
    rnonlinearSpatialMatrix%Dalpha(2,1) = 1.0_DP
    rnonlinearSpatialMatrix%Dalpha(1,2) = 1.0_DP
    rnonlinearSpatialMatrix%Dalpha(2,2) = 1.0_DP
    rnonlinearSpatialMatrix%Dnewton(2,1) = real(1-rphysicsPrimal%iequation,DP)
      
  end subroutine

  ! ***************************************************************************
  
!<subroutine>

  subroutine cc_spaceTimeMatVec (rproblem, rspaceTimeMatrix, rx, rd, cx, cy, &
      cfilter, dnorm, bprintRes)

!<description>
  ! This routine performs a matrix-vector multiplication with the
  ! system matrix A defined by rspaceTimeMatrix.
  !    rd  :=  cx A(p_rsolution) rx  +  cy rd
  ! If rspaceTimeMatrix does not specify an evaluation point for th nonlinearity
  ! in A(.), the routine calculates
  !    rd  :=  cx A(rx) rx  +  cy rd
  ! rd is overwritten by the result.
!</description>

!<input>
  ! A problem structure that provides information about matrices on all
  ! levels as well as temporary vectors.
  type(t_problem), intent(INOUT), target :: rproblem

  ! A t_ccoptSpaceTimeMatrix structure defining the space-time matrix.
  type(t_ccoptSpaceTimeMatrix), intent(IN) :: rspaceTimeMatrix

  ! A space-time vector defining the current solution.
  type(t_spacetimeVector), intent(IN) :: rx
  
  ! Type of filter to apply to the vectors. A combination of SPTID_FILTER_xxxx
  ! flags that specifies which type of filter is to be applied to the
  ! output vector during the matriux vector multiplication.
  integer(I32), intent(IN) :: cfilter
  
  ! If set to TRUE and dnorm is present, too, the residuals are printed to the
  ! terminal.
  logical, intent(IN), optional :: bprintRes
!</input>

!<inputoutput>
  ! A second space-time vector.
  type(t_spacetimeVector), intent(INOUT) :: rd
  
  ! Multiplication factor for rx
  real(DP), intent(IN) :: cx
  
  ! Multiplication factor for rd
  real(DP), intent(IN) :: cy
!</inputoutput>

!<output>
  ! OPTIONAL: If specified, returns the $l_2$-norm of rd.
  real(DP), intent(OUT), optional :: dnorm
!<output>

!</subroutine>

    ! local variables
    integer :: ieqTime,ilevel,icp
    type(t_vectorBlock) :: rtempVectorD
    type(t_vectorBlock), dimension(3) :: rtempVector
    type(t_vectorBlock), dimension(3) :: rtempVectorEval
    type(t_blockDiscretisation), pointer :: p_rdiscr
    real(DP) :: dtheta,dnormpart
    type(t_matrixBlock) :: rblockTemp
    type(t_nonlinearSpatialMatrix) :: rnonlinearSpatialMatrix
    type(t_ccoptSpaceTimeDiscretisation), pointer :: p_rspaceTimeDiscr
    real(DP) :: dtstep,dtime
    type(t_discreteBC) :: rdiscreteBC
    type(t_discreteFBC) :: rdiscreteFBC
    
    ! DEBUG!!!
    real(DP), dimension(:), pointer :: p_Dx1,p_Dx2,p_Dx3,p_Db
    real(DP), dimension(:), pointer :: p_DxE1,p_DxE2,p_DxE3
    
    ! Pointer to the space time discretisation
    p_rspaceTimeDiscr => rspaceTimeMatrix%p_rspaceTimeDiscr
    
    ! Level of the discretisation
    ilevel = p_rspaceTimeDiscr%ilevel

    ! Theta-scheme identifier.
    ! =1: implicit Euler.
    ! =0.5: Crank Nicolson
    dtheta = rproblem%rtimedependence%dtimeStepTheta
   
    ! Initialise the boundary conditions
    call bcasm_initDiscreteBC(rdiscreteBC)
    call bcasm_initDiscreteFBC(rdiscreteFBC)

    ! Create a temp vector that contains the part of rd which is to be modified.
    p_rdiscr => p_rspaceTimeDiscr%p_RlevelInfo%rdiscretisation
    call lsysbl_createVecBlockByDiscr (p_rdiscr,rtempVectorD,.false.)
    
    ! The vector will be a defect vector. Assign the boundary conditions so
    ! that we can implement them.
    call lsysbl_assignDiscreteBC(rtempVectorD,rdiscreteBC)
    call lsysbl_assignDiscreteFBC(rtempVectorD,rdiscreteFBC)
    
    ! Create a temp vector for the X-vectors at timestep i-1, i and i+1.
    call lsysbl_createVecBlockByDiscr (p_rdiscr,rtempVector(1),.false.)
    call lsysbl_createVecBlockByDiscr (p_rdiscr,rtempVector(2),.false.)
    call lsysbl_createVecBlockByDiscr (p_rdiscr,rtempVector(3),.false.)
    
    ! Get the parts of the X-vector which are to be modified at first --
    ! subvector 1, 2 and 3. rtempVector(1) contains the 'previous' solution,
    ! rtempVector(2) the 'current' and rtempVector(3) the 'next' one.
    call sptivec_getTimestepData(rx, 1+0, rtempVector(2))
    
    ! If necesary, multiply the rtempVectorX. We have to take a -1 into
    ! account as the actual matrix multiplication routine cc_assembleDefect
    ! introduces another -1!
    if (cx .ne. -1.0_DP) then
      call lsysbl_scaleVector (rtempVector(2),-cx)
    end if
      
    ! Now what's the evaluation point where to evaluate the nonlinearity/ies?
    ! If the structure defines an evaluation point, we take that one, so
    ! we need another temp vector that holds the evaluation point.
    ! If not, we take the vector rx. This is archived by creating new
    ! vectors that share their information with rx.
    if (associated(rspaceTimeMatrix%p_rsolution)) then
      call lsysbl_createVecBlockByDiscr (p_rdiscr,rtempVectorEval(1),.false.)
      call lsysbl_createVecBlockByDiscr (p_rdiscr,rtempVectorEval(2),.false.)
      call lsysbl_createVecBlockByDiscr (p_rdiscr,rtempVectorEval(3),.false.)
    else
      call lsysbl_duplicateVector (rtempVector(1),rtempVectorEval(1),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      call lsysbl_duplicateVector (rtempVector(2),rtempVectorEval(2),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      call lsysbl_duplicateVector (rtempVector(3),rtempVectorEval(3),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    end if

    ! If a nonlinearity is involved, rtempVectorEval(1) contains the 'previous',
    ! rtempVectorEval(2) the 'current' and rtempVectorEval(3) the 'next'
    ! solution (relative to the current time step) where to evaluate
    ! the nonlinearity.
    if (associated(rspaceTimeMatrix%p_rsolution)) then
      call sptivec_getTimestepData(rspaceTimeMatrix%p_rsolution, &
          1+0, rtempVectorEval(2))
    end if
    
    ! DEBUG!!!
    call lsysbl_getbase_double (rtempVector(1),p_Dx1)
    call lsysbl_getbase_double (rtempVector(2),p_Dx2)
    call lsysbl_getbase_double (rtempVector(3),p_Dx3)
    call lsysbl_getbase_double (rtempVectorD,p_Db)
    call lsysbl_getbase_double (rtempVectorEval(1),p_DxE1)
    call lsysbl_getbase_double (rtempVectorEval(2),p_DxE2)
    call lsysbl_getbase_double (rtempVectorEval(3),p_DxE3)
    
    ! If dnorm is specified, clear it.
    if (present(dnorm)) then
      dnorm = 0.0_DP
    end if

    dtstep = p_rspaceTimeDiscr%rtimeDiscr%dtstep
    
    ! Loop through the substeps
    
    do ieqTime = 0,p_rspaceTimeDiscr%NEQtime-1
      
      ! Depending on the time step scheme, initialise the current time.
      select case (p_rspaceTimeDiscr%rtimeDiscr%ctype)
      case (TDISCR_ONESTEPTHETA)
        dtime = rproblem%rtimedependence%dtimeInit + ieqTime*dtstep
      case (TDISCR_DG0)
        dtime = rproblem%rtimedependence%dtimeInit + (real(ieqTime,DP)+0.5_DP)*dtstep
      end select
      
      ! Get the part of rd which is to be modified.
      if (cy .ne. 0.0_DP) then
        call sptivec_getTimestepData(rd, 1+ieqTime, rtempVectorD)
        
        ! If cy <> 1, multiply rtempVectorD by that.
        if (cy .ne. 1.0_DP) then
          call lsysbl_scaleVector (rtempVectorD,cy)
        end if
      else
        call lsysbl_clearVector (rtempVectorD)
      end if
      
      if (ieqTime .ne. p_rspaceTimeDiscr%NEQtime-1) then
      
        ! Read the solution of the 'next' timestep and the 'next' evaluation
        ! point.

        call sptivec_getTimestepData(rx, 1+ieqTime+1, rtempVector(3))

        ! If necesary, multiply the rtempVectorX. We have to take a -1 into
        ! account as the actual matrix multiplication routine cc_assembleDefect
        ! introduces another -1!
        if (cx .ne. -1.0_DP) then
          call lsysbl_scaleVector (rtempVector(3),-cx)
        end if

        if (associated(rspaceTimeMatrix%p_rsolution)) then
          call sptivec_getTimestepData(rspaceTimeMatrix%p_rsolution, &
              1+ieqTime+1, rtempVectorEval(3))

        end if

      end if

      ! The first and last substep is a little bit special concerning
      ! the matrix!
      if (ieqTime .eq. 0) then
        
        ! We are in the first substep. Here, we have to handle the following
        ! part of the supersystem:
        !
        !  ( A11 A12   0   0 ... )  ( x1 )  =  ( f1 )
        !  ( ... ... ... ... ... )  ( .. )     ( .. )
        !
        ! So we have to compute:
        !
        !  d1  :=  f1  -  A11 x1  -  A12 x2
        !
        ! -----
        
        ! The diagonal matrix.
      
        ! Set up the matrix weights of that submatrix.
        call cc_initNonlinMatrix (rnonlinearSpatialMatrix,rproblem,&
            rspaceTimeMatrix%p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
            rspaceTimeMatrix%p_rspaceTimeDiscr%p_rlevelInfo%rstaticInfo)
        call cc_setupMatrixWeights (rproblem,rspaceTimeMatrix,dtheta,&
          ieqTime,0,rnonlinearSpatialMatrix)
          
        ! Subtract: rd = rd - A11 x1
        call cc_assembleDefect (rnonlinearSpatialMatrix,rtempVector(2),rtempVectorD,&
            1.0_DP,rtempVectorEval(1),rtempVectorEval(2),rtempVectorEval(3))

        ! -----
      
        ! Create the matrix
        !   A12  :=  -M + dt*dtheta*[-nu\Laplace u + u \grad u]
        ! and subtract A12 x2 from rd.

        ! Set up the matrix weights of that submatrix.
        call cc_initNonlinMatrix (rnonlinearSpatialMatrix,rproblem,&
            rspaceTimeMatrix%p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
            rspaceTimeMatrix%p_rspaceTimeDiscr%p_rlevelInfo%rstaticInfo)
        call cc_setupMatrixWeights (rproblem,rspaceTimeMatrix,dtheta,&
          ieqTime,1,rnonlinearSpatialMatrix)

        ! Subtract: rd = rd - A12 x2
        call cc_assembleDefect (rnonlinearSpatialMatrix,rtempVector(3),rtempVectorD,&
            1.0_DP,rtempVectorEval(1),rtempVectorEval(2),rtempVectorEval(3))

        ! Release the block mass matrix.
        call lsysbl_releaseMatrix (rblockTemp)

      else if (ieqTime .lt. p_rspaceTimeDiscr%NEQtime-1) then

        ! We are sonewhere in the middle of the matrix. There is a substep
        ! ieqTime+1 and a substep ieqTime-1!  Here, we have to handle the following
        ! part of the supersystem:
        !
        !  ( ... ...   ... ...   ... )  ( .. )     ( .. )
        !  ( ... Aii-1 Aii Aii+1 ... )  ( xi )  =  ( fi )
        !  ( ... ...   ... ...   ... )  ( .. )     ( .. )
        !
        ! So we have to compute:
        !
        !  dn  :=  fn  -  Aii-1 xi-1  -  Aii xi  -  Aii+1 xi+1
        
        ! -----
        
        ! Create the matrix
        !   Aii-1 := -M + dt*dtheta*[-nu\Laplace u + u \grad u]
        ! and include that into the global matrix for the primal velocity.

        ! Set up the matrix weights of that submatrix.
        call cc_initNonlinMatrix (rnonlinearSpatialMatrix,rproblem,&
            rspaceTimeMatrix%p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
            rspaceTimeMatrix%p_rspaceTimeDiscr%p_rlevelInfo%rstaticInfo)
        call cc_setupMatrixWeights (rproblem,rspaceTimeMatrix,dtheta,&
          ieqTime,-1,rnonlinearSpatialMatrix)
            
        ! Subtract: rd = rd - Aii-1 xi-1.
        ! Note that at this point, the nonlinearity must be evaluated
        ! at xi due to the discretisation scheme!!!
        call cc_assembleDefect (rnonlinearSpatialMatrix,rtempVector(1),rtempVectorD,&
            1.0_DP,rtempVectorEval(1),rtempVectorEval(2),rtempVectorEval(3))

        ! Release the block mass matrix.
        call lsysbl_releaseMatrix (rblockTemp)

        ! -----

        ! Now the diagonal matrix.
      
        ! Assemble the nonlinear defect.
      
        ! Set up the matrix weights of that submatrix.
        call cc_initNonlinMatrix (rnonlinearSpatialMatrix,rproblem,&
            rspaceTimeMatrix%p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
            rspaceTimeMatrix%p_rspaceTimeDiscr%p_rlevelInfo%rstaticInfo)
        call cc_setupMatrixWeights (rproblem,rspaceTimeMatrix,dtheta,&
          ieqTime,0,rnonlinearSpatialMatrix)

        ! Subtract: rd = rd - Aii xi
        call cc_assembleDefect (rnonlinearSpatialMatrix,rtempVector(2),rtempVectorD,&
            1.0_DP,rtempVectorEval(1),rtempVectorEval(2),rtempVectorEval(3))
            
        ! -----
        
        ! Create the matrix
        !   Aii+1 := -M + dt*dtheta*[-nu\Laplace u + u \grad u]
        ! and include that into the global matrix for the dual velocity.

        ! Set up the matrix weights of that submatrix.
        call cc_initNonlinMatrix (rnonlinearSpatialMatrix,rproblem,&
            rspaceTimeMatrix%p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
            rspaceTimeMatrix%p_rspaceTimeDiscr%p_rlevelInfo%rstaticInfo)
        call cc_setupMatrixWeights (rproblem,rspaceTimeMatrix,dtheta,&
          ieqTime,1,rnonlinearSpatialMatrix)
          
        ! Subtract: rd = rd - Aii+1 xi+1
        call cc_assembleDefect (rnonlinearSpatialMatrix,rtempVector(3),rtempVectorD,&
            1.0_DP,rtempVectorEval(1),rtempVectorEval(2),rtempVectorEval(3))
        
        ! Release the block mass matrix.
        call lsysbl_releaseMatrix (rblockTemp)

      else
        
        ! We are in the last substep. Here, we have to handle the following
        ! part of the supersystem:
        !
        !  ( ... ... ... ...   ... )  ( .. )     ( .. )
        !  ( ... ...   0 Ann-1 Ann )  ( xn )  =  ( fn )
        !
        ! So we have to compute:
        !
        !  dn  :=  fn  -  Ann-1 xn-1  -  Ann xn
        
        ! -----
        
        ! Create the matrix
        !   Ann-1 = -M + dt*dtheta*[-nu\Laplace u + u \grad u]
        ! and include that into the global matrix for the dual velocity.

        ! Set up the matrix weights of that submatrix.
        call cc_initNonlinMatrix (rnonlinearSpatialMatrix,rproblem,&
            rspaceTimeMatrix%p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
            rspaceTimeMatrix%p_rspaceTimeDiscr%p_rlevelInfo%rstaticInfo)
        call cc_setupMatrixWeights (rproblem,rspaceTimeMatrix,dtheta,&
          ieqTime,-1,rnonlinearSpatialMatrix)
          
        ! Subtract: rd = rd - Ann-1 xn-1
        ! Note that at this point, the nonlinearity must be evaluated
        ! at xn due to the discretisation scheme!!!
        call cc_assembleDefect (rnonlinearSpatialMatrix,rtempVector(1),rtempVectorD,&
            1.0_DP,rtempVectorEval(1),rtempVectorEval(2),rtempVectorEval(3))
     
        ! -----
        
        ! Now the diagonal matrix.
      
        ! Assemble the nonlinear defect.
      
        ! Set up the matrix weights of that submatrix.
        call cc_initNonlinMatrix (rnonlinearSpatialMatrix,rproblem,&
            rspaceTimeMatrix%p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
            rspaceTimeMatrix%p_rspaceTimeDiscr%p_rlevelInfo%rstaticInfo)
        call cc_setupMatrixWeights (rproblem,rspaceTimeMatrix,dtheta,&
          ieqTime,0,rnonlinearSpatialMatrix)

        ! Subtract: rd = rd - Ann xn
        call cc_assembleDefect (rnonlinearSpatialMatrix,rtempVector(2),rtempVectorD,&
            1.0_DP,rtempVectorEval(1),rtempVectorEval(2),rtempVectorEval(3))
      
      end if

      ! Implement the BC`s?
      if (iand(cfilter,SPTID_FILTER_BCDEF) .ne. 0) then

        ! Discretise the boundary conditions at the new point in time.
        call bcasm_clearDiscreteBC(rdiscreteBC)
        call bcasm_clearDiscreteFBC(rdiscreteFBC)
        call cc_assembleBDconditions (rproblem,dtime,&
            p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
            CCDISCBC_PRIMALDUAL,rdiscreteBC,rproblem%rcollection)
        call cc_assembleFBDconditions (rproblem,dtime,&
            p_rspaceTimeDiscr%p_rlevelInfo%rdiscretisation,&
            CCDISCBC_PRIMALDUAL,rdiscreteFBC,rproblem%rcollection)

        ! Implement the boundary conditions into the defect.
        call vecfil_discreteBCdef (rtempVectorD)
        call vecfil_discreteFBCdef (rtempVectorD)
      end if

      if ((ieqTime .eq. 0) .and. (iand(cfilter,SPTID_FILTER_ICDEF) .ne. 0)) then
        ! Implement the initial conditions into the defect.
        call tbc_implementInitCondDefSingle (p_rspaceTimeDiscr, rtempVectorD)
      end if

      if ((ieqTime .eq. p_rspaceTimeDiscr%NEQtime-1) .and. &
          (iand(cfilter,SPTID_FILTER_TCDEF) .ne. 0)) then
        ! Implement the initial conditions into the defect.
        call tbc_implementTermCondDefSingle (p_rspaceTimeDiscr, rtempVectorD)
      end if
      
      ! Save the defect vector back to rd.
      call sptivec_setTimestepData(rd, 1+ieqTime, rtempVectorD)
      
      ! If dnorm is specified, calculate the norm of the sub-defect vector and
      ! add it to dnorm.
      if (present(dnorm)) then
        dnormpart = lsysbl_vectorNorm(rtempVectorD,LINALG_NORML2)**2
        dnorm = dnorm + dnormpart
        
        if (present(bprintRes)) then
          if (bprintRes) then
            call output_line ('||D_'//trim(sys_siL(1+ieqTime,10))//'|| = '//&
                trim(sys_sdEL(sqrt(dnormpart),10)) )
            do icp=1,rtempVectorD%nblocks
              call output_line ('  ||D_'//&
                  trim(sys_siL(1+ieqTime,10))//'^'//trim(sys_siL(icp,2))&
                  //'|| = '//&
                  trim(sys_sdEL(&
                      lsyssc_vectorNorm(rtempVectorD%RvectorBlock(icp),LINALG_NORML2),10)) )
            end do
          end if
        end if
      end if
      
      ! Cycle the solution vectors and the evaluation points: 1 <- 2 <- 3.
      call lsysbl_copyVector (rtempVector(2),rtempVector(1))
      call lsysbl_copyVector (rtempVector(3),rtempVector(2))

      if (associated(rspaceTimeMatrix%p_rsolution)) then
        call lsysbl_copyVector (rtempVectorEval(2),rtempVectorEval(1))
        call lsysbl_copyVector (rtempVectorEval(3),rtempVectorEval(2))
      end if
      
      ! Now, the 3rd subvector is again free to take the next vector.
      
    end do
    
    ! If dnorm is specified, normalise it.
    ! It was calculated from rspaceTimeDiscr%niterations+1 subvectors.
    if (present(dnorm)) then
      dnorm = sqrt(dnorm / real(p_rspaceTimeDiscr%NEQtime,DP))
    end if
    
    ! Release the temp vectors.
    call lsysbl_releaseVector (rtempVectorEval(3))
    call lsysbl_releaseVector (rtempVectorEval(2))
    call lsysbl_releaseVector (rtempVectorEval(1))
    call lsysbl_releaseVector (rtempVector(3))
    call lsysbl_releaseVector (rtempVector(2))
    call lsysbl_releaseVector (rtempVector(1))
    call lsysbl_releaseVector (rtempVectorD)
    
    ! Release the BC's again.
    call bcasm_releaseDiscreteFBC(rdiscreteFBC)
    call bcasm_releaseDiscreteBC(rdiscreteBC)
    
  end subroutine
   
  ! ***************************************************************************
  
!<subroutine>

  subroutine stlin_disableSubmatrix (rnonlinearSpatialMatrix,irow,icolumn)

!<description>
  ! Disables a subbklock in the nonlinear matrix rnonlinearSpatialMatrix.
  ! All weights of the correspopnding subblock are set to 0.
!</description>

!<input>
  ! The row/column of the submatrix to be disabled.
  integer :: irow,icolumn
!</input>

!<inputoutput>
  ! A t_nonlinearSpatialMatrix structure that defines the shape of the core
  ! equation. The weights that specify the submatrices of a small 6x6
  ! block matrix system are initialised depending on the position
  ! specified by isubstep and nsubsteps.
  !
  ! The structure must have been initialised with smva_initNonlinMatrix!
  type(t_nonlinearSpatialMatrix), intent(inout) :: rnonlinearSpatialMatrix
!</inputoutput>

!</subroutine>

    ! Clear the coefficients
    rnonlinearSpatialMatrix%Diota(irow,icolumn) = 0.0_DP
    rnonlinearSpatialMatrix%Dalpha(irow,icolumn) = 0.0_DP
    rnonlinearSpatialMatrix%Dtheta(irow,icolumn) = 0.0_DP
    rnonlinearSpatialMatrix%Dgamma(irow,icolumn) = 0.0_DP
    rnonlinearSpatialMatrix%Dnewton(irow,icolumn) = 0.0_DP
    rnonlinearSpatialMatrix%DgammaT(irow,icolumn) = 0.0_DP
    rnonlinearSpatialMatrix%Dnewton2(irow,icolumn) = 0.0_DP
    rnonlinearSpatialMatrix%DgammaT2(irow,icolumn) = 0.0_DP
    rnonlinearSpatialMatrix%DnewtonT(irow,icolumn) = 0.0_DP
    rnonlinearSpatialMatrix%Deta(irow,icolumn) = 0.0_DP
    rnonlinearSpatialMatrix%Dtau(irow,icolumn) = 0.0_DP
    rnonlinearSpatialMatrix%Dkappa(irow,icolumn) = 0.0_DP

  end subroutine


  ! ***************************************************************************
  
!<subroutine>

  subroutine cc_generateInitCondRHS (rproblem,rspaceTimeDiscr,rx,rb)

!<description>
  ! Generates the RHS vector used for the initial condition.
!</description>

!<input>
  ! Problem structure of the main problem.
  type(t_problem), intent(INOUT) :: rproblem
  
  ! Space-time discretisation structure of the level of the RHS.
  type(t_ccoptSpaceTimeDiscretisation), intent(IN) :: rspaceTimeDiscr
  
  ! Solution vector containing the solution of the first timestep.
  type(t_vectorBlock), intent(in) :: rx
!</input>

!<inputoutput>
  ! Vector that receives the RHS for the initial condition.
  ! The vector must have been initialised. The vector data is overwritten
  ! by the RHS data.
  type(t_vectorBlock), intent(inout) :: rb
!</inputoutput>

!</subroutine>

    ! local variables
    type(t_nonlinearSpatialMatrix) :: rnonlinearSpatialMatrix
    type(t_problem_lvl), pointer :: p_rlevelInfo
    logical :: bconvectionExplicit
    real(DP) :: dtstep, dequationType
    real(dp), dimension(:), pointer :: p_Ddata
    
!    ! If the following constant is set from 1.0 to 0.0, the primal system is
!    ! decoupled from the dual system!
!    real(DP), parameter :: dprimalDualCoupling = 1.0_DP
!
!    ! If the following constant is set from 1.0 to 0.0, the dual system is
!    ! decoupled from the primal system!
!    real(DP), parameter :: ddualPrimalCoupling = 1.0_DP
!
!    ! If the following parameter is set from 1.0 to 0.0, the time coupling
!    ! is disabled, resulting in a stationary simulation in every timestep.
!    real(DP), parameter :: dtimeCoupling = 1.0_DP

    real(DP) :: dprimalDualCoupling,ddualPrimalCoupling,dtimeCoupling
    
    call parlst_getvalue_double (rproblem%rparamList,'DEBUG',&
        'dprimalDualCoupling',dprimalDualCoupling,1.0_DP)

    call parlst_getvalue_double (rproblem%rparamList,'DEBUG',&
        'ddualPrimalCoupling',ddualPrimalCoupling,1.0_DP)

    call parlst_getvalue_double (rproblem%rparamList,'DEBUG',&
        'dtimeCoupling',dtimeCoupling,1.0_DP)

    ! The initial condition is implemented as:
    !
    !   (M/dt + A) y_0  =  b_0 := (M/dt + A) y^0
    !
    ! i.e. we take the solution vector of the 0th timestep. multiply it by the
    ! (Navier--)Stokes equation and what we receive is the RHS for the
    ! terminal condition. We only have to take care of bondary conditions.
    
    call lsysbl_getbase_double (rb,p_Ddata)

    ! Set up the basic components of the Navier--Stokes matrix
    p_rlevelInfo => rspaceTimeDiscr%p_rlevelInfo

    call cc_initNonlinMatrix (rnonlinearSpatialMatrix,rproblem,&
        p_rlevelInfo%rdiscretisation,p_rlevelInfo%rstaticInfo)

    ! Disable the submatrices for the dual solution and the coupling.
    ! We only want to generate the RHS for the primal solution.
    call stlin_disableSubmatrix (rnonlinearSpatialMatrix,2,1)
    call stlin_disableSubmatrix (rnonlinearSpatialMatrix,1,2)
    call stlin_disableSubmatrix (rnonlinearSpatialMatrix,2,2)

    ! Change the sign of dupsam2 for a consistent stabilisation.
    ! Reason: The stablisation is added to the dual operator by the SD/
    ! EOJ stabilisation in the following way:
    !
    !    ... - (u grad lamda + dupsam2*stabilisation) + ... = rhs
    !
    ! We want to *add* the stabilisation, so we have to introduce a "-" sign
    ! in dupsam2 to get
    !
    !    ... - (u grad lamda) - (-dupsam2*stabilisation) + ... = rhs
    ! <=>
    !    ... - (u grad lamda) + dupsam2*stabilisation + ... = rhs
    
    rnonlinearSpatialMatrix%dupsam2 = -rnonlinearSpatialMatrix%dupsam2
    
    ! Set up the matrix weights the matrix in the 0th timestep.
    ! We set up only the stuff for the primal equation; for setting up
    ! the RHS, there is no dual equation and also no coupling between the
    ! primal and dual solutions.
    bconvectionExplicit = rproblem%roptcontrol%iconvectionExplicit .ne. 0

    dtstep = rspaceTimeDiscr%rtimeDiscr%dtstep
    
    dequationType = 1.0_DP
    if (rproblem%roptcontrol%ispaceTimeFormulation .ne. 0) &
      dequationType = -1.0_DP
      
    rnonlinearSpatialMatrix%Dalpha(1,1) = dtimeCoupling * 1.0_DP/dtstep
    rnonlinearSpatialMatrix%Dtheta(1,1) = 1.0_DP
    
    if (.not. bconvectionExplicit) then
      rnonlinearSpatialMatrix%Dgamma(1,1) = real(1-rproblem%rphysicsPrimal%iequation,DP)
    end if

    rnonlinearSpatialMatrix%Deta(1,1) = 1.0_DP
    rnonlinearSpatialMatrix%Dtau(1,1) = 1.0_DP
        
    ! Create by substraction: rd = 0*rd - (- A11 x1) = A11 x1
    ! Clear the primal RHS for that purpose.
    call lsyssc_clearVector(rb%RvectorBlock(1))
    call lsyssc_clearVector(rb%RvectorBlock(2))
    call lsyssc_clearVector(rb%RvectorBlock(3))
    call cc_assembleDefect (rnonlinearSpatialMatrix,rx,rb,-1.0_DP,rx,rx,rx)

  end subroutine

!  ! ***************************************************************************
!
!!<subroutine>
!
!  SUBROUTINE cc_assembleSpaceTimeRHS (rproblem, rspaceTimeDiscr, rb, &
!      rtempvector1, rtempvector2, rtempvector3, bimplementBC)
!
!!<description>
!  ! Assembles the space-time RHS vector rb. Bondary conditions are NOT
!  ! implemented!
!  !
!  ! Note: rproblem%rtimedependence%dtime will be undefined at the end of
!  ! this routine!
!!</description>
!
!!<input>
!  ! A problem structure that provides information on all
!  ! levels as well as temporary vectors.
!  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem
!
!  ! A t_ccoptSpaceTimeDiscretisation structure defining the discretisation of the
!  ! coupled space-time system.
!  TYPE(t_ccoptSpaceTimeDiscretisation), INTENT(IN) :: rspaceTimeDiscr
!!</input>
!
!!<inputoutput>
!  ! A temporary vector in the size of a spatial vector.
!  TYPE(t_vectorBlock), INTENT(INOUT) :: rtempVector1
!
!  ! A second temporary vector in the size of a spatial vector.
!  TYPE(t_vectorBlock), INTENT(INOUT) :: rtempVector2
!
!  ! A third temporary vector in the size of a spatial vector.
!  TYPE(t_vectorBlock), INTENT(INOUT) :: rtempVector3
!
!  ! A space-time vector that receives the RHS.
!  TYPE(t_spacetimeVector), INTENT(INOUT) :: rb
!
!  ! Whether to implement boundary conditions into the RHS or not.
!  LOGICAL, INTENT(IN) :: bimplementBC
!!</inputoutput>
!
!!</subroutine>
!
!    ! local variables
!    INTEGER :: isubstep
!    REAL(DP) :: dtheta
!    TYPE(t_nonlinearSpatialMatrix) :: rnonlinearSpatialMatrix
!
!    ! A temporary vector for the creation of the RHS.
!    TYPE(t_vectorBlock) :: rtempVectorRHS
!
!    REAL(DP), DIMENSION(:),POINTER :: p_Dx, p_Db, p_Dd, p_Drhs
!
!    ! Theta-scheme identifier.
!    ! =1: impliciz Euler.
!    ! =0.5: Crank Nicolson
!    dtheta = rproblem%rtimedependence%dtimeStepTheta
!
!    ! ----------------------------------------------------------------------
!    ! Generate the global RHS vector
!
!    CALL lsysbl_getbase_double (rtempVector1,p_Dx)
!    CALL lsysbl_getbase_double (rtempVector2,p_Db)
!    CALL lsysbl_getbase_double (rtempVector3,p_Dd)
!
!    ! Assemble 1st RHS vector in X temp vector.
!    CALL generateRHS (rproblem,0,rb%ntimesteps,rspaceTimeDiscr%rtimeDiscr%dtstep,&
!        rtempVector1, .TRUE., .FALSE.)
!
!    ! Assemble the 2nd RHS vector in the RHS temp vector
!    CALL generateRHS (rproblem,1,rb%ntimesteps,rspaceTimeDiscr%rtimeDiscr%dtstep,&
!        rtempVector2, .TRUE., .FALSE.)
!
!    ! Assemble the 3rd RHS vector in the defect temp vector
!    IF (rspaceTimeDiscr%rtimeDiscr%nintervals .GE. 2) THEN
!      CALL generateRHS (rproblem,2,rb%ntimesteps,rspaceTimeDiscr%rtimeDiscr%dtstep,&
!          rtempVector3, .TRUE., .FALSE.)
!    ELSE
!      CALL lsysbl_copyVector (rtempVector2,rtempVector3)
!    END IF
!
!    ! Create a copy of the X temp vector (RHS0). That vector will be
!    ! our destination vector for assembling the RHS in all timesteps.
!    CALL lsysbl_copyVector (rtempVector1,rtempVectorRHS)
!
!    ! DEBUG!!!
!    CALL lsysbl_getbase_double (rtempVectorRHS,p_Drhs)
!
!    ! RHS 0,1,2 -> 1-2-3
!
!    DO isubstep = 0,rspaceTimeDiscr%rtimeDiscr%nintervals
!
!      IF (isubstep .EQ. 0) THEN
!
!        ! Primal RHS comes from rtempVector1. The dual from the
!        ! isubstep+1'th RHS in rtempVector2.
!        !
!        ! primal RHS(0) = PRIMALRHS(0)
!        ! dual RHS(0)   = THETA*DUALRHS(0) + (1-THETA)*DUALRHS(1)
!
!        CALL lsyssc_copyVector (rtempVector1%RvectorBlock(1),rtempVectorRHS%RvectorBlock(1))
!        CALL lsyssc_copyVector (rtempVector1%RvectorBlock(2),rtempVectorRHS%RvectorBlock(2))
!        CALL lsyssc_copyVector (rtempVector1%RvectorBlock(3),rtempVectorRHS%RvectorBlock(3))
!
!        CALL lsyssc_vectorLinearComb (&
!            rtempVector1%RvectorBlock(4),rtempVector2%RvectorBlock(4),&
!            dtheta,(1.0_DP-dtheta),&
!            rtempVectorRHS%RvectorBlock(4))
!        CALL lsyssc_vectorLinearComb (&
!            rtempVector1%RvectorBlock(5),rtempVector2%RvectorBlock(5),&
!            dtheta,(1.0_DP-dtheta),&
!            rtempVectorRHS%RvectorBlock(5))
!        ! Pressure is fully implicit
!        CALL lsyssc_vectorLinearComb (&
!            rtempVector1%RvectorBlock(6),rtempVector2%RvectorBlock(6),&
!            dtheta,(1.0_DP-dtheta),&
!            rtempVectorRHS%RvectorBlock(6))
!
!        ! In the 0'th timestep, there is no RHS in the dual equation!
!        ! That is because of the initial condition, which fixes the primal solution
!        ! => dual solution has no influence on the primal one
!        ! => setting up a dual RHS in not meaningful as the dual RHS cannot
!        !    influence the primal solution
!        !CALL lsyssc_clearVector (rtempVectorRHS%RvectorBlock(4))
!        !CALL lsyssc_clearVector (rtempVectorRHS%RvectorBlock(5))
!        !CALL lsyssc_clearVector (rtempVectorRHS%RvectorBlock(6))
!
!      ELSE IF (isubstep .LT. rspaceTimeDiscr%rtimeDiscr%nintervals) THEN
!
!        ! We are somewhere 'in the middle'.
!        !
!        ! Dual RHS comes from rtempVector3. The primal from the
!        ! isubstep-1'th RHS.
!        !
!        ! primal RHS(0) = THETA*PRIMALRHS(0) + (1-THETA)*PRIMALRHS(-1)
!        ! dual RHS(0)   = THETA*DUALRHS(0) + (1-THETA)*DUALRHS(1)
!
!        CALL lsyssc_vectorLinearComb (&
!            rtempVector1%RvectorBlock(1),rtempVector2%RvectorBlock(1),&
!            (1.0_DP-dtheta),dtheta,&
!            rtempVectorRHS%RvectorBlock(1))
!        CALL lsyssc_vectorLinearComb (&
!            rtempVector1%RvectorBlock(2),rtempVector2%RvectorBlock(2),&
!            (1.0_DP-dtheta),dtheta,&
!            rtempVectorRHS%RvectorBlock(2))
!        ! Pressure is fully implicit
!        CALL lsyssc_vectorLinearComb (&
!            rtempVector1%RvectorBlock(3),rtempVector2%RvectorBlock(3),&
!            (1.0_DP-dtheta),dtheta,&
!            rtempVectorRHS%RvectorBlock(3))
!
!        CALL lsyssc_vectorLinearComb (&
!            rtempVector2%RvectorBlock(4),rtempVector3%RvectorBlock(4),&
!            dtheta,(1.0_DP-dtheta),&
!            rtempVectorRHS%RvectorBlock(4))
!        CALL lsyssc_vectorLinearComb (&
!            rtempVector2%RvectorBlock(5),rtempVector3%RvectorBlock(5),&
!            dtheta,(1.0_DP-dtheta),&
!            rtempVectorRHS%RvectorBlock(5))
!        ! Pressure is fully implicit
!        CALL lsyssc_vectorLinearComb (&
!            rtempVector2%RvectorBlock(6),rtempVector3%RvectorBlock(6),&
!            dtheta,(1.0_DP-dtheta),&
!            rtempVectorRHS%RvectorBlock(6))
!
!        IF (isubstep .LT. rspaceTimeDiscr%rtimeDiscr%nintervals-1) THEN
!          ! Shift the RHS vectors and generate the RHS for the next time step.
!          ! (Yes, I know, this could probably be solved more elegant without copying anything
!          ! using a ring buffer ^^)
!          CALL lsysbl_copyVector(rtempVector2,rtempVector1)
!          CALL lsysbl_copyVector(rtempVector3,rtempVector2)
!          CALL generateRHS (rproblem,isubstep+2,&
!              rspaceTimeDiscr%rtimeDiscr%nintervals,&
!              rspaceTimeDiscr%rtimeDiscr%dtstep,rtempVector3, .TRUE., .FALSE.)
!        END IF
!
!      ELSE
!
!        ! We are 'at the end'.
!        !
!        ! Dual RHS comes from rtempVector3. The primal from the
!        ! isubstep-1'th RHS and rtempVector3.
!        !
!        ! primal RHS(0) = THETA*PRIMALRHS(0) + (1-THETA)*PRIMALRHS(-1)
!        ! dual RHS(0)   = DUALRHS(0)
!
!        CALL lsyssc_vectorLinearComb (&
!            rtempVector2%RvectorBlock(1),rtempVector3%RvectorBlock(1),&
!            (1.0_DP-dtheta),dtheta,&
!            rtempVectorRHS%RvectorBlock(1))
!        CALL lsyssc_vectorLinearComb (&
!            rtempVector2%RvectorBlock(2),rtempVector3%RvectorBlock(2),&
!            (1.0_DP-dtheta),dtheta,&
!            rtempVectorRHS%RvectorBlock(2))
!        ! Pressure is fully implicit
!        CALL lsyssc_vectorLinearComb (&
!            rtempVector2%RvectorBlock(3),rtempVector3%RvectorBlock(3),&
!            (1.0_DP-dtheta),dtheta,&
!            rtempVectorRHS%RvectorBlock(3))
!
!        !CALL generateRHS (rproblem,isubstep+1,rspaceTimeDiscr%niterations,&
!        !    rtempVector3, .TRUE., .FALSE.)
!
!        CALL lsyssc_copyVector (rtempVector3%RvectorBlock(4),rtempVectorRHS%RvectorBlock(4))
!        CALL lsyssc_copyVector (rtempVector3%RvectorBlock(5),rtempVectorRHS%RvectorBlock(5))
!        CALL lsyssc_copyVector (rtempVector3%RvectorBlock(6),rtempVectorRHS%RvectorBlock(6))
!
!        ! Multiply the last RHS of the dual equation -z by gamma, that's it.
!        CALL lsyssc_scaleVector (rtempVectorRHS%RvectorBlock(4),rspaceTimeDiscr%dgammaC)
!        CALL lsyssc_scaleVector (rtempVectorRHS%RvectorBlock(5),rspaceTimeDiscr%dgammaC)
!
!      END IF
!
!      ! Implement the boundary conditions into the RHS vector
!      IF (bimplementBC) THEN
!        CALL generateRHS (rproblem,isubstep,&
!            rspaceTimeDiscr%rtimeDiscr%nintervals,&
!            rspaceTimeDiscr%rtimeDiscr%dtstep, rtempVectorRHS, .FALSE., .TRUE.)
!      END IF
!
!      ! Save the RHS.
!      CALL sptivec_setTimestepData(rb, isubstep, rtempVectorRHS)
!
!    END DO
!
!    ! Release the temp vector for generating the RHS.
!    CALL lsysbl_releaseVector (rtempVectorRHS)
!
!  CONTAINS
!
!    SUBROUTINE generateRHS (rproblem,isubstep,nsubsteps,dtstep,&
!        rvector, bgenerate, bincludeBC)
!
!    ! Generate the RHS vector of timestep isubstep and/or include boundary
!    ! conditions.
!
!    ! Problem structure.
!    TYPE(t_problem), INTENT(INOUT) :: rproblem
!
!    ! Number of the substep where to generate the RHS vector
!    INTEGER, INTENT(IN) :: isubstep
!
!    ! Total number of substeps
!    INTEGER, INTENT(IN) :: nsubsteps
!
!    ! Length od one timestep
!    REAL(DP), INTENT(IN) :: dtstep
!
!    ! Destination vector
!    TYPE(t_vectorBlock), INTENT(INOUT) :: rvector
!
!    ! Whether to generate the RHS vector or not
!    LOGICAL, INTENT(IN) :: bgenerate
!
!    ! Whether to include boundary conditions
!    LOGICAL, INTENT(IN) :: bincludeBC
!
!    ! DEBUG!!!
!    REAL(DP), DIMENSION(:), POINTER :: p_Ddata
!
!      ! DEBUG!!!
!      CALL lsysbl_getbase_double (rvector,p_Ddata)
!
!      ! Set the time where we are at the moment
!      rproblem%rtimedependence%dtime = &
!          rproblem%rtimedependence%dtimeInit + isubstep*dtstep
!      rproblem%rtimedependence%itimestep = isubstep
!
!      ! Assemble the RHS?
!      IF (bgenerate) THEN
!
!        ! Generate the RHS of that point in time into the vector.
!        CALL cc_generateBasicRHS (rproblem,rvector)
!
!      END IF
!
!      ! Include BC's?
!      IF (bincludeBC) THEN
!
!        ! Initialise the collection for the assembly process with callback routines.
!        ! Basically, this stores the simulation time in the collection if the
!        ! simulation is nonstationary.
!        CALL cc_initCollectForAssembly (rproblem,rproblem%rcollection)
!
!        ! Discretise the boundary conditions at the new point in time --
!        ! if the boundary conditions are nonconstant in time!
!        IF (collct_getvalue_int (rproblem%rcollection,'IBOUNDARY') .NE. 0) THEN
!          CALL cc_updateDiscreteBC (rproblem)
!        END IF
!
!        ! Implement the boundary conditions into the RHS.
!        ! This is done *after* multiplying -z by GAMMA or dtstep, resp.,
!        ! as Dirichlet values mustn't be multiplied with GAMMA!
!        CALL vecfil_discreteBCsol (rvector)
!        CALL vecfil_discreteFBCsol (rvector)
!
!        ! Clean up the collection (as we are done with the assembly, that's it.
!        CALL cc_doneCollectForAssembly (rproblem,rproblem%rcollection)
!
!      END IF
!
!    END SUBROUTINE
!
!  END SUBROUTINE

!betrachte die GMV's. Entkoppelt gibt's 2. Ordnung in der Geschwindigkeit
!(Druck=0), Gekoppelt nicht :(
  
end module
