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
!# Auxiliary routines:
!#
!# 1.) cc_setupMatrixWeights
!#     -> Initialise the matrix weights of a submatrix in the global space-time
!#        matrix.
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
!#
!#                   Newton iteration:
!#                   A = A(y_i) = -nu*Laplace(.) + y_i*grad(.) + grad(y_i)*(.)
!#                   N = N(y_i) = -nu*Laplace(.) - y_i*grad(.) + grad(y_i)*(.)
!#                   R = R(l_i) = -M + l_i*grad(.) - grad(l_i)*(.)
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
!#  [-M/dt ]                       | [M/dt + A] [-B] [ 1/a M   ]      |
!#                                 |                                  |
!#                                 | [-B^t    ]                       |
!#                                 |                                  |
!#                                 | [ R      ]      [M/dt + N ] [-B] |                  [-M/dt   ]
!#                                 |   ^                              |
!#                                 | from the        [-B^t     ]      |
!#                                 | RHS                              |
!#  -------------------------------+----------------------------------+---------------------------------
!#                                 | [-M/dt   ]                       |  [M/dt + A] [-B] [ 1/a M  ]      
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
!#  [-M/dt+(1-T)A(y_0)]                                   | [M/dt + T A(y_1) ] [-B ] [ 1/a M           ]       |                                            
!#                                                        |                                                    |                                            
!#                                                        | [-B^t            ]                                 |                                            
!#                                                        |                                                    |                                            
!#                                                        | [-T M            ]       [ M/dt + T N(y_1) ] [-B ] |                          [-M/dt+(1-T)N(y_1)]
!#                                                        |   ^                                                |                                        ^    
!#                                                        | from the                 [ -B^t            ]       |                                        |    
!#                                                        | RHS                                                |                                    !correct!    
!#  ------------------------------------------------------+----------------------------------------------------+---------------------------------------------------
!#                                                        | [-M/dt+(1-T)A(y_1)]                                |  [M/dt + T A(y_2)] [-B ] [ 1/a M           ]       
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

MODULE spacetimelinearsystem

  USE fsystem
  USE storage
  USE linearsolver
  USE boundary
  USE bilinearformevaluation
  USE linearformevaluation
  USE cubature
  USE matrixfilters
  USE vectorfilters
  USE bcassembly
  USE triangulation
  USE spatialdiscretisation
  USE coarsegridcorrection
  USE spdiscprojection
  USE nonlinearsolver
  USE paramlist
  USE linearsolverautoinitialise
  USE matrixrestriction
  USE paramlist
  USE timestepping
  USE l2projection
  
  USE collection
  USE convection
    
  USE cc2dmediumm2basic
  USE cc2dmedium_callback

  USE cc2dmediumm2nonlinearcore
  USE cc2dmediumm2nonlinearcoreinit
  USE cc2dmediumm2timeanalysis
  USE cc2dmediumm2boundary
  USE cc2dmediumm2discretisation
  USE cc2dmediumm2matvecassembly
  
  USE timediscretisation
  USE spacetimediscretisation
  USE spacetimevectors
  USE timeboundaryconditions
  
  USE matrixio
    
  IMPLICIT NONE

!<types>

!<typeblock>

  ! Defines the basic shape of the supersystem which realises the coupling
  ! between all timesteps.
  TYPE t_ccoptSpaceTimeMatrix
  
    ! Type of the matrix.
    ! =0: Standard system matrix.
    ! =1: Newton matrix
    INTEGER :: cmatrixType = 0
  
    ! Pointer to the space time discretisation that corresponds to this matrix.
    TYPE(t_ccoptSpaceTimeDiscretisation), POINTER :: p_rspaceTimeDiscretisation => NULL()
  
    ! Pointer to a space-time solution vector that defines the point
    ! where the nonlinearity is evaluated when applying or calculating
    ! matrices.
    ! If this points to NULL(), there is no global solution specified,
    ! so the matrix assembly/application routines use the vector that
    ! is specified in their input parameters.
    TYPE(t_spacetimeVector), POINTER :: p_rsolution => NULL()
    
  END TYPE

!</typeblock>

!</types>

!<constants>

!<constantblock description="Filter-ID's of the filters to be applied in a matrix vector multiplication">

  ! No filter
  INTEGER(I32), PARAMETER :: SPTID_FILTER_NONE  = 0

  ! Filter the output vector for boundary conditions in the defect.
  INTEGER(I32), PARAMETER :: SPTID_FILTER_BCDEF = 2**0
  
  ! Filter the output vector for initial conditions in the defect
  INTEGER(I32), PARAMETER :: SPTID_FILTER_ICDEF = 2**1

  ! Filter the output vector for terminal conditions in the defect
  INTEGER(I32), PARAMETER :: SPTID_FILTER_TCDEF = 2**2
  
  ! Apply all filters that are typically applied to a defect vector
  INTEGER(I32), PARAMETER :: SPTID_FILTER_DEFECT = SPTID_FILTER_BCDEF + &
                                                   SPTID_FILTER_ICDEF + &
                                                   SPTID_FILTER_TCDEF

!</constantblock>

!</constants>

CONTAINS

  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE cc_setupMatrixWeights (rproblem,rspaceTimeMatrix,dtheta,&
      isubstep,irelpos,rmatrixComponents)

!<description>
  ! This routine sets up the matrix weights in the rmatrixComponents structure
  ! according to the position of the corresponding submatrix.
  ! isubstep defines the timestep which is to be tackled -- which corresponds
  ! to the row in the supermatrix.
  ! irelpos specifies the column in the supermatrix relative to the diagonal --
  ! thus =0 means the diagonal, =1 the submatrix above the diagonal and =-1
  ! the submatrix below the diagonal.
  !
  ! The matrix weights in rmatrixComponents are initialised based on this
  ! information such that the 'assembleMatrix' and 'assembleDefect' routines of 
  ! the core equation will build the system matrix / defect at that position 
  ! in the supermatrix.
  !
  ! The routine does not initialise the pointers to the basic matrices/discretisation
  ! structures in the rmatrixComponents structure. This has to be done by
  ! the caller!
!</description>

!<input>
  ! Problem structure
  TYPE(t_problem), INTENT(IN) :: rproblem
  
  ! A t_ccoptSpaceTimeDiscretisation structure defining the discretisation of the
  ! coupled space-time matrix.
  TYPE(t_ccoptSpaceTimeMatrix), INTENT(IN), TARGET :: rspaceTimeMatrix

  ! Theta scheme identifier.
  ! = 0.5: Crank-Nicolson.
  ! = 1.0: Explicit Euler
  REAL(DP), INTENT(IN) :: dtheta
  
  ! Substep in the time-dependent simulation = row in the supermatrix.
  ! Range 0..nsubsteps
  INTEGER, INTENT(IN) :: isubstep
  
  ! Specifies the column in the supermatrix relative to the diagonal.
  ! =0: set matrix weights for the diagonal.
  INTEGER, INTENT(IN) :: irelpos
!</input>

!<inputoutput>
  ! A t_ccmatrixComponents structure that defines the shape of the core
  ! equation. The weights that specify the submatrices of a small 6x6 
  ! block matrix system are initialised depending on the position
  ! specified by isubstep and nsubsteps.
  TYPE(t_ccmatrixComponents), INTENT(INOUT) :: rmatrixComponents
!</inputoutput>

!</subroutine>

    ! If the following constant is set from 1.0 to 0.0, the primal system is
    ! decoupled from the dual system!
    REAL(DP), PARAMETER :: dprimalDualCoupling = 1.0_DP
    
    ! If the following constant is set from 1.0 to 0.0, the dual system is
    ! decoupled from the primal system!
    REAL(DP), PARAMETER :: ddualPrimalCoupling = 1.0_DP
    
    ! If the following parameter is set from 1.0 to 0.0, the terminal
    ! condition between the primal and dual equation is decoupled, i.e.
    ! the dual equation gets independent from the primal one.
    REAL(DP), PARAMETER :: dterminalCondDecoupled = 1.0_DP
    
    ! If the following parameter is set from 1.0 to 0.0, the time coupling
    ! is disabled, resulting in a stationary simulation in every timestep.
    REAL(DP), PARAMETER :: dtimeCoupling = 1.0_DP
    
    ! This constant defines the type of equation. There are two equivalent
    ! formulations of the dual equation which only differs in the sign
    ! of the dual velocity.
    ! A constant of "1" here means to use the formulation with "y-z" in
    ! the RHS of the dual equation, while a constant of "-1" means to use the
    ! formulation with "-(y-z)" there.
    ! Note: If this is changed, a "-" sign must be implemented / removed from
    ! the RHS, too!
    REAL(DP) :: dequationType
    
    ! Pointer to the space time discretisation structure.
    TYPE(t_ccoptSpaceTimeDiscretisation), POINTER :: p_rspaceTimeDiscr
    REAL(DP) :: dnewton
    REAL(DP) :: dsmooth
    REAL(DP) :: dtstep
    LOGICAL :: bconvectionExplicit
    
    p_rspaceTimeDiscr => rspaceTimeMatrix%p_rspaceTimeDiscretisation
    
    dequationType = 1.0_DP
    IF (rproblem%roptcontrol%ispaceTimeFormulation .NE. 0) &
      dequationType = -1.0_DP

    ! Treat the convection explicitely?
    bconvectionExplicit = rproblem%roptcontrol%iconvectionExplicit .NE. 0

    ! What's the matrix type we should set up? If we have to set up a
    ! Newton matrix, we put dnewton to 1.0 so the Newton part in the
    ! primal velocity is assembled.
    dnewton = 0.0_DP
    IF (rspaceTimeMatrix%cmatrixType .EQ. 1) THEN
      IF (rproblem%iequation .EQ. 0) THEN
        ! Newton is only to be assembled in Navier-Stokes!
        dnewton = 1.0_DP
      END IF
      !dsmooth = 1.0_DP
    ELSE 
      ! Activate smoothing in the last timestep.
      !dsmooth = 1.0_DP
    END IF
    
    ! Note: Introducing the Navier-Stokes operator to smooth the last time step
    ! is a difficile task. The best combination found up to now is:
    ! Always switch on the 'smoothing' (dsmooth=1) of the terminal condition 
    ! to get a proper matrix and filter the defect of the last timestep to 
    ! impose the terminal condition!
    !dsmooth = 1.0_DP
    !
    ! We activate the smoothing in the matrix contrary to the filter.
    ! If the filter is on (thus weakly imposing the terminal conditions
    ! to the defect by filtering), we activate the Navier-Stokes-Operator in 
    ! the last timestep. If the 'smoothing' is on (thus the defect vectors in 
    ! the last timestep of the dual solution are not forced to be zero), we
    ! deactivate the Navier-Stokes-operator.

    ! Probably activate smoothing in the last timestep.
    dsmooth = REAL(1-p_rspaceTimeDiscr%itypeTerminalCondition,DP)

    dtstep = p_rspaceTimeDiscr%rtimeDiscr%dtstep

    ! The first and last substep is a little bit special concerning
    ! the matrix!
    IF (isubstep .EQ. 0) THEN
      
      ! We are in the first substep
      !
      ! Specify which vectors should be used for evaluating the nonlinearity.
      ! Aii-1 = undefined.
      ! Aii = Aii(x_i,lambda_i).
      ! Aii+1 = Aii+1(x_i,lambda_i+1)
    
      IF (irelpos .EQ. 0) THEN
      
        ! The diagonal matrix.
        rmatrixComponents%iprimalSol = 2
        rmatrixComponents%idualSol = 2

        rmatrixComponents%diota1 = 1.0_DP
        rmatrixComponents%diota2 = 0.0_DP

        rmatrixComponents%dkappa1 = 1.0_DP
        rmatrixComponents%dkappa2 = 0.0_DP
        
        rmatrixComponents%dalpha1 = 0.0_DP
        rmatrixComponents%dalpha2 = dtimeCoupling * 1.0_DP/dtstep
        
        rmatrixComponents%dtheta1 = 0.0_DP
        rmatrixComponents%dtheta2 = dtheta
        
        IF (.NOT. bconvectionExplicit) THEN
        
          rmatrixComponents%dgamma1 = 0.0_DP
          rmatrixComponents%dgamma2 = &
              - dtheta * REAL(1-rproblem%iequation,DP)
          
          rmatrixComponents%dnewton1 = 0.0_DP
          rmatrixComponents%dnewton2 = &
                dtheta * REAL(1-rproblem%iequation,DP)
                
        ELSE

          rmatrixComponents%dgamma1 = 0.0_DP
          rmatrixComponents%dgamma2 = 0.0_DP
          
          rmatrixComponents%dnewton1 = 0.0_DP
          rmatrixComponents%dnewton2 = 0.0_DP

        END IF

        rmatrixComponents%deta1 = 0.0_DP
        rmatrixComponents%deta2 = 1.0_DP
        
        rmatrixComponents%dtau1 = 0.0_DP
        rmatrixComponents%dtau2 = 1.0_DP
        
        rmatrixComponents%dmu1 = 0.0_DP
        rmatrixComponents%dmu2 = ddualPrimalCoupling * &
            dtheta * (-dequationType)

        IF (.NOT. bconvectionExplicit) THEN
        
          ! In the 0'th timestep, there is no RHS in the dual equation
          ! and therefore no coupling between the primal and dual solution!
          ! That is because of the initial condition, which fixes the primal solution
          ! => dual solution has no influence on the primal one
          ! Therefore, the following weights must be commented out, otherwise
          ! the solver cannot converge in the 0'th timestep!
          ! Instead, the dual solution of this 0'th timestep only depends
          ! on the dual solution of the 1st timestep and is computed assuming that
          ! the initial condition meets the target flow (y_0 = z_0).
          !
          ! rmatrixComponents%dmu1 = dprimalDualCoupling * &
          !     dtheta * 1.0_DP / p_rspaceTimeDiscr%dalphaC

          IF (dnewton .EQ. 0.0_DP) THEN
            rmatrixComponents%dr21 = 0.0_DP
            rmatrixComponents%dr22 = 0.0_DP
          ELSE
            rmatrixComponents%dr21 = ddualPrimalCoupling * &
                dtheta * ( dequationType)
            rmatrixComponents%dr22 = ddualPrimalCoupling * &
                dtheta * (-dequationType)
          END IF
          
        ELSE
        
          rmatrixComponents%dr21 = 0.0_DP
          rmatrixComponents%dr22 = 0.0_DP

        END IF
          
                    
      ELSE IF (irelpos .EQ. 1) THEN
      
        ! Offdiagonal matrix on the right of the diagonal.
        !
        ! Note that at this point, the nonlinearity must be evaluated
        ! at xi due to the discretisation scheme!!!
        !
        ! WARNING: For a very strange situation, taking xi here (which is said
        ! to be the correct evaluation point from the theory) does not lead
        ! to quadratic convergence in the Newton. Taking xi+1 does!?!
        ! So we take xi+1, although the theory tells us to take xi!

        rmatrixComponents%iprimalSol = 2
        rmatrixComponents%idualSol = 3

        ! Create the matrix
        !   -M + dt*dtheta*[-nu\Laplace u + u \grad u]

        rmatrixComponents%diota1 = 0.0_DP
        rmatrixComponents%diota2 = 0.0_DP

        rmatrixComponents%dkappa1 = 0.0_DP
        rmatrixComponents%dkappa2 = 0.0_DP
        
        rmatrixComponents%dalpha1 = 0.0_DP
        rmatrixComponents%dalpha2 = dtimeCoupling * (-1.0_DP)/dtstep
        
        rmatrixComponents%dtheta1 = 0.0_DP
        rmatrixComponents%dtheta2 = (1.0_DP-dtheta) 
        
        IF (.NOT. bconvectionExplicit) THEN
        
          rmatrixComponents%dgamma1 = 0.0_DP
          rmatrixComponents%dgamma2 = &
              - (1.0_DP-dtheta) * REAL(1-rproblem%iequation,DP)
          
          rmatrixComponents%dnewton1 = 0.0_DP
          rmatrixComponents%dnewton2 = &
                (1.0_DP-dtheta) * REAL(1-rproblem%iequation,DP)
        ELSE

          rmatrixComponents%dgamma1 = 0.0_DP
          rmatrixComponents%dgamma2 = &
              - dtheta * REAL(1-rproblem%iequation,DP)
          
          rmatrixComponents%dnewton1 = 0.0_DP
          rmatrixComponents%dnewton2 = &
                dtheta * REAL(1-rproblem%iequation,DP)

        END IF

        rmatrixComponents%deta1 = 0.0_DP
        rmatrixComponents%deta2 = 0.0_DP
        
        rmatrixComponents%dtau1 = 0.0_DP
        rmatrixComponents%dtau2 = 0.0_DP
        
        rmatrixComponents%dmu1 = 0.0_DP
        rmatrixComponents%dmu2 = ddualPrimalCoupling * &
            (-dequationType) * (1.0_DP-dtheta)
            
        IF (.NOT. bconvectionExplicit) THEN

          IF (dnewton .EQ. 0.0_DP) THEN
            rmatrixComponents%dr21 = 0.0_DP
            rmatrixComponents%dr22 = 0.0_DP
          ELSE
            rmatrixComponents%dr21 = ddualPrimalCoupling * &
                (1.0_DP-dtheta) * ( dequationType)
            rmatrixComponents%dr22 = ddualPrimalCoupling * &
                (1.0_DP-dtheta) * (-dequationType)
          END IF
        
        ELSE

          rmatrixComponents%dr21 = 0.0_DP
          rmatrixComponents%dr22 = 0.0_DP
        
        END IF
            
      END IF
    
    ELSE IF (isubstep .LT. p_rspaceTimeDiscr%NEQtime-1) THEN
      
      ! We are sonewhere in the middle of the matrix. There is a substep
      ! isubstep+1 and a substep isubstep-1!
      
      ! -----
      
      ! Specify which vectors should be used for evaluating the nonlinearity.
      ! Aii-1 = Aii+1(x_i-1,lambda_i-1).
      ! Aii = Aii(x_i,lambda_i).
      ! Aii+1 = Aii+1(x_i,lambda_i+1)

      IF (irelpos .EQ. -1) THEN
      
        ! Matrix on the left of the diagonal.

        rmatrixComponents%iprimalSol = 1
        rmatrixComponents%idualSol = 1

        ! Create the matrix
        !   -M + dt*dtheta*[-nu\Laplace u + u \grad u]

        rmatrixComponents%diota1 = 0.0_DP
        rmatrixComponents%diota2 = 0.0_DP

        rmatrixComponents%dkappa1 = 0.0_DP
        rmatrixComponents%dkappa2 = 0.0_DP
        
        rmatrixComponents%dalpha1 = dtimeCoupling * (-1.0_DP)/dtstep
        rmatrixComponents%dalpha2 = 0.0_DP
        
        rmatrixComponents%dtheta1 = (1.0_DP-dtheta) 
        rmatrixComponents%dtheta2 = 0.0_DP
        
        IF (.NOT. bconvectionExplicit) THEN

          rmatrixComponents%dgamma1 = &
              (1.0_DP-dtheta) * REAL(1-rproblem%iequation,DP)
          rmatrixComponents%dgamma2 = 0.0_DP
          
          rmatrixComponents%dnewton1 = &
              (1.0_DP-dtheta) * dnewton
          rmatrixComponents%dnewton2 = 0.0_DP

        ELSE
        
          rmatrixComponents%dgamma1 = &
              dtheta * REAL(1-rproblem%iequation,DP)
          rmatrixComponents%dgamma2 = 0.0_DP
          
          rmatrixComponents%dnewton1 = &
              dtheta * dnewton
          rmatrixComponents%dnewton2 = 0.0_DP
        
        END IF

        rmatrixComponents%deta1 = 0.0_DP
        rmatrixComponents%deta2 = 0.0_DP
        
        rmatrixComponents%dtau1 = 0.0_DP
        rmatrixComponents%dtau2 = 0.0_DP
        
        rmatrixComponents%dmu1 = dprimalDualCoupling * &
            dequationType * (1.0_DP-dtheta) / p_rspaceTimeDiscr%dalphaC
        rmatrixComponents%dmu2 = 0.0_DP

        rmatrixComponents%dr21 = 0.0_DP
        rmatrixComponents%dr22 = 0.0_DP

      ELSE IF (irelpos .EQ. 0) THEN    

        ! The diagonal matrix.
        IF (.NOT. bconvectionExplicit) THEN
          rmatrixComponents%iprimalSol = 2
          rmatrixComponents%idualSol = 2
        ELSE
          rmatrixComponents%iprimalSol = 2
          rmatrixComponents%idualSol = 3
        END IF

        rmatrixComponents%diota1 = 0.0_DP
        rmatrixComponents%diota2 = 0.0_DP

        rmatrixComponents%dkappa1 = 0.0_DP
        rmatrixComponents%dkappa2 = 0.0_DP
        
        rmatrixComponents%dalpha1 = dtimeCoupling * 1.0_DP/dtstep
        rmatrixComponents%dalpha2 = dtimeCoupling * 1.0_DP/dtstep
        
        rmatrixComponents%dtheta1 = dtheta
        rmatrixComponents%dtheta2 = dtheta
        
        IF (.NOT. bconvectionExplicit) THEN

          rmatrixComponents%dgamma1 = &
              dtheta * REAL(1-rproblem%iequation,DP)
          rmatrixComponents%dgamma2 = &
              - dtheta * REAL(1-rproblem%iequation,DP)
          
          rmatrixComponents%dnewton1 = dtheta * dnewton
          rmatrixComponents%dnewton2 = &
                dtheta * REAL(1-rproblem%iequation,DP)
                
        ELSE
        
          rmatrixComponents%dgamma1 = 0.0_DP
          rmatrixComponents%dgamma2 = 0.0_DP
          
          rmatrixComponents%dnewton1 = 0.0_DP
          rmatrixComponents%dnewton2 = 0.0_DP
        
        END IF

        rmatrixComponents%deta1 = 1.0_DP
        rmatrixComponents%deta2 = 1.0_DP
        
        rmatrixComponents%dtau1 = 1.0_DP
        rmatrixComponents%dtau2 = 1.0_DP
        
        rmatrixComponents%dmu1 = dprimalDualCoupling * &
            dequationType * dtheta * 1.0_DP / p_rspaceTimeDiscr%dalphaC
        rmatrixComponents%dmu2 = ddualPrimalCoupling * &
            (-dequationType) * dtheta 
            
        IF (.NOT. bconvectionExplicit) THEN

          IF (dnewton .EQ. 0.0_DP) THEN
            rmatrixComponents%dr21 = 0.0_DP
            rmatrixComponents%dr22 = 0.0_DP
          ELSE
            rmatrixComponents%dr21 = ddualPrimalCoupling * &
                ( dequationType) * dtheta 
            rmatrixComponents%dr22 = ddualPrimalCoupling * &
                (-dequationType) * dtheta 
          END IF
          
        ELSE
        
          IF (dnewton .EQ. 0.0_DP) THEN
            rmatrixComponents%dr21 = 0.0_DP
            rmatrixComponents%dr22 = 0.0_DP
          ELSE
            rmatrixComponents%dr21 = ddualPrimalCoupling * &
                ( dequationType) * dtheta 
            rmatrixComponents%dr22 = ddualPrimalCoupling * &
                (-dequationType) * dtheta 
          END IF
        
        END IF

      ELSE IF (irelpos .EQ. 1) THEN
            
        ! Matrix on the right of the diagonal.
        !
        ! Note that at this point, the nonlinearity must be evaluated
        ! at xi due to the discretisation scheme!!!
        !
        ! WARNING: For a very strange situation, taking xi here (which is said
        ! to be the correct evaluation point from the theory) does not lead
        ! to quadratic convergence in the Newton. Taking xi+1 does!?!
        ! So we take xi+1, although the theory tells us to take xi!
        rmatrixComponents%iprimalSol = 2
        rmatrixComponents%idualSol = 3

        ! Create the matrix
        !   -M + dt*dtheta*[-nu\Laplace u + u \grad u]
        rmatrixComponents%diota1 = 0.0_DP
        rmatrixComponents%diota2 = 0.0_DP

        rmatrixComponents%dkappa1 = 0.0_DP
        rmatrixComponents%dkappa2 = 0.0_DP
        
        rmatrixComponents%dalpha1 = 0.0_DP
        rmatrixComponents%dalpha2 = dtimeCoupling * (-1.0_DP)/dtstep
        
        rmatrixComponents%dtheta1 = 0.0_DP
        rmatrixComponents%dtheta2 = (1.0_DP-dtheta) 
        
        IF (.NOT. bconvectionExplicit) THEN

          rmatrixComponents%dgamma1 = 0.0_DP
          rmatrixComponents%dgamma2 = &
              - (1.0_DP-dtheta) * REAL(1-rproblem%iequation,DP)
          
          rmatrixComponents%dnewton1 = 0.0_DP
          rmatrixComponents%dnewton2 = &
              (1.0_DP-dtheta) * REAL(1-rproblem%iequation,DP)
              
        ELSE
        
          rmatrixComponents%dgamma1 = 0.0_DP
          rmatrixComponents%dgamma2 = &
              - dtheta * REAL(1-rproblem%iequation,DP)
          
          rmatrixComponents%dnewton1 = 0.0_DP
          rmatrixComponents%dnewton2 = &
                dtheta * REAL(1-rproblem%iequation,DP)
        
        END IF

        rmatrixComponents%deta1 = 0.0_DP
        rmatrixComponents%deta2 = 0.0_DP
        
        rmatrixComponents%dtau1 = 0.0_DP
        rmatrixComponents%dtau2 = 0.0_DP
        
        rmatrixComponents%dmu1 = 0.0_DP
        rmatrixComponents%dmu2 = ddualPrimalCoupling * &
            (-dequationType) * (1.0_DP-dtheta) 
            
        IF (.NOT. bconvectionExplicit) THEN
        
          IF (dnewton .EQ. 0.0_DP) THEN
            rmatrixComponents%dr21 = 0.0_DP
            rmatrixComponents%dr22 = 0.0_DP
          ELSE
            rmatrixComponents%dr21 = ddualPrimalCoupling * &
                (1.0_DP-dtheta) * ( dequationType)
            rmatrixComponents%dr22 = ddualPrimalCoupling * &
                (1.0_DP-dtheta) * (-dequationType)
          END IF
          
        ELSE
        
          rmatrixComponents%dr21 = 0.0_DP
          rmatrixComponents%dr22 = 0.0_DP
          
        END IF
            
      END IF
    
    ELSE
    
      ! We are in the last substep
      
      ! Specify which vectors should be used for evaluating the nonlinearity.
      ! Aii-1 = Aii+1(x_i-1,lambda_i-1).
      ! Aii = Aii(x_i,lambda_i).
      ! Aii+1 = Aii+1(x_i,lambda_i+1)

      IF (irelpos .EQ. -1) THEN
      
        ! Matrix on the left of the diagonal.
        
        rmatrixComponents%iprimalSol = 1
        rmatrixComponents%idualSol = 1

        ! Create the matrix
        !   -M + dt*dtheta*[-nu\Laplace u + u \grad u]

        rmatrixComponents%diota1 = 0.0_DP
        rmatrixComponents%diota2 = 0.0_DP

        rmatrixComponents%dkappa1 = 0.0_DP
        rmatrixComponents%dkappa2 = 0.0_DP
        
        rmatrixComponents%dalpha1 = dtimeCoupling * (-1.0_DP)/dtstep
        rmatrixComponents%dalpha2 = 0.0_DP
        
        rmatrixComponents%dtheta1 = (1.0_DP-dtheta) 
        rmatrixComponents%dtheta2 = 0.0_DP
        
        IF (.NOT. bconvectionExplicit) THEN

          rmatrixComponents%dgamma1 = &
              (1.0_DP-dtheta) * REAL(1-rproblem%iequation,DP)
          rmatrixComponents%dgamma2 = 0.0_DP
          
          rmatrixComponents%dnewton1 = (1.0_DP-dtheta) * dnewton
          rmatrixComponents%dnewton2 = 0.0_DP
          
        ELSE
        
          rmatrixComponents%dgamma1 = &
              dtheta * REAL(1-rproblem%iequation,DP)
          rmatrixComponents%dgamma2 = 0.0_DP
          
          rmatrixComponents%dnewton1 = dtheta * dnewton
          rmatrixComponents%dnewton2 = 0.0_DP
          
        END IF

        rmatrixComponents%deta1 = 0.0_DP
        rmatrixComponents%deta2 = 0.0_DP
        
        rmatrixComponents%dtau1 = 0.0_DP
        rmatrixComponents%dtau2 = 0.0_DP
        
        rmatrixComponents%dmu1 = dprimalDualCoupling * &
            dequationType * (1.0_DP-dtheta) / p_rspaceTimeDiscr%dalphaC
        rmatrixComponents%dmu2 = 0.0_DP

        rmatrixComponents%dr21 = 0.0_DP
        rmatrixComponents%dr22 = 0.0_DP

      ELSE IF (irelpos .EQ. 0) THEN    

        ! The diagonal matrix.
        rmatrixComponents%iprimalSol = 2
        rmatrixComponents%idualSol = 2

        rmatrixComponents%diota1 = 0.0_DP
        rmatrixComponents%diota2 = 0.0_DP ! (1.0_DP-dsmooth) !0.0_DP

        rmatrixComponents%dkappa1 = 0.0_DP
        rmatrixComponents%dkappa2 = (1.0_DP-dsmooth) !0.0_DP
        
        ! Current formulation:
        ! -gamma*M*y + (M+dt*nu*L)*lambda = -gamma*z
        
        rmatrixComponents%dalpha1 = dtimeCoupling * 1.0_DP/dtstep
        rmatrixComponents%dalpha2 = 1.0_DP
        !(1.0_DP-dsmooth) + dsmooth ! * 1.0_DP/dtstep
        ! No 'time coupling' here; because of the terminal condition,
        ! the mass matrix resembles not the time dependence!
        
        rmatrixComponents%dtheta1 = dtheta
        rmatrixComponents%dtheta2 = dsmooth * dtheta * dtstep
        
        IF (.NOT. bconvectionExplicit) THEN

          rmatrixComponents%dgamma1 = &
              dtheta * REAL(1-rproblem%iequation,DP)
          rmatrixComponents%dgamma2 = &
              - dsmooth * dtheta * REAL(1-rproblem%iequation,DP) * dtstep
          
          rmatrixComponents%dnewton1 = dtheta * dnewton
          rmatrixComponents%dnewton2 = &
                dsmooth * dtheta * REAL(1-rproblem%iequation,DP) * dtstep
                
        ELSE
        
          rmatrixComponents%dgamma1 = 0.0_DP
          rmatrixComponents%dgamma2 = 0.0_DP
          
          rmatrixComponents%dnewton1 = 0.0_DP
          rmatrixComponents%dnewton2 = 0.0_DP
        
        END IF        

        rmatrixComponents%deta1 = 1.0_DP
        rmatrixComponents%deta2 = dsmooth * dtstep
        
        rmatrixComponents%dtau1 = 1.0_DP
        rmatrixComponents%dtau2 = dsmooth
        
        rmatrixComponents%dmu1 = dprimalDualCoupling * &
            dequationType * dtheta * 1.0_DP / p_rspaceTimeDiscr%dalphaC
        rmatrixComponents%dmu2 = ddualPrimalCoupling * &
            (-dequationType) * dtheta * p_rspaceTimeDiscr%dgammaC
            
        IF (.NOT. bconvectionExplicit) THEN
        
          ! Weight the mass matrix by GAMMA instead of delta(T).
          ! That's the only difference to the implementation above!
          IF (dnewton .EQ. 0.0_DP) THEN
            rmatrixComponents%dr21 = 0.0_DP
            rmatrixComponents%dr22 = 0.0_DP
          ELSE
            rmatrixComponents%dr21 = ddualPrimalCoupling * &
                ( dequationType) * dtheta * dtstep
            rmatrixComponents%dr22 = ddualPrimalCoupling * &
                (-dequationType) * dtheta * dtstep
          END IF

        ELSE
        
          rmatrixComponents%dr21 = 0.0_DP
          rmatrixComponents%dr22 = 0.0_DP

        END IF

      END IF        
        
!      IF (irelpos .EQ. -1) THEN
!      
!        ! Matrix on the left of the diagonal.
!        !
!        ! Create the matrix
!        !   -M + dt*dtheta*[-nu\Laplace u + u \grad u]
!
!        rmatrixComponents%diota1 = 0.0_DP
!        rmatrixComponents%diota2 = 0.0_DP
!
!        rmatrixComponents%dkappa1 = 0.0_DP
!        rmatrixComponents%dkappa2 = 0.0_DP
!        
!        rmatrixComponents%dalpha1 = -1.0_DP
!        rmatrixComponents%dalpha2 = 0.0_DP
!        
!        rmatrixComponents%dtheta1 = (1.0_DP-dtheta) * p_rspaceTimeDiscr%dtstep
!        rmatrixComponents%dtheta2 = 0.0_DP
!        
!        rmatrixComponents%dgamma1 = &
!            (1.0_DP-dtheta) * p_rspaceTimeDiscr%dtstep * REAL(1-rproblem%iequation,DP)
!        rmatrixComponents%dgamma2 = 0.0_DP
!        
!        rmatrixComponents%dnewton1 = 0.0_DP
!        rmatrixComponents%dnewton2 = 0.0_DP
!
!        rmatrixComponents%deta1 = 0.0_DP
!        rmatrixComponents%deta2 = 0.0_DP
!        
!        rmatrixComponents%dtau1 = 0.0_DP
!        rmatrixComponents%dtau2 = 0.0_DP
!        
!        rmatrixComponents%dmu1 = dprimalDualCoupling * &
!            p_rspaceTimeDiscr%dtstep * (1.0_DP-dtheta) / p_rspaceTimeDiscr%dalphaC
!        rmatrixComponents%dmu2 = 0.0_DP
!
!      ELSE IF (irelpos .EQ. 0) THEN    
!
!        ! The diagonal matrix.
!
!        rmatrixComponents%diota1 = 0.0_DP
!        rmatrixComponents%diota2 = 0.0_DP
!
!        rmatrixComponents%dkappa1 = 0.0_DP
!        rmatrixComponents%dkappa2 = 0.0_DP
!        
!        rmatrixComponents%dalpha1 = 1.0_DP
!        rmatrixComponents%dalpha2 = 1.0_DP
!        
!        rmatrixComponents%dtheta1 = dtheta * p_rspaceTimeDiscr%dtstep
!        rmatrixComponents%dtheta2 = dtheta * p_rspaceTimeDiscr%dtstep
!        
!        rmatrixComponents%dgamma1 = &
!            dtheta * p_rspaceTimeDiscr%dtstep * REAL(1-rproblem%iequation,DP)
!        rmatrixComponents%dgamma2 = &
!            - dtheta * p_rspaceTimeDiscr%dtstep * REAL(1-rproblem%iequation,DP)
!        
!        rmatrixComponents%dnewton1 = 0.0_DP
!        rmatrixComponents%dnewton2 = &
!              dtheta * p_rspaceTimeDiscr%dtstep * REAL(1-rproblem%iequation,DP)
!
!        rmatrixComponents%deta1 = p_rspaceTimeDiscr%dtstep
!        rmatrixComponents%deta2 = p_rspaceTimeDiscr%dtstep
!        
!        rmatrixComponents%dtau1 = 1.0_DP
!        rmatrixComponents%dtau2 = 1.0_DP
!        
!        rmatrixComponents%dmu1 = dprimalDualCoupling * &
!            dtheta * p_rspaceTimeDiscr%dtstep / p_rspaceTimeDiscr%dalphaC
!        rmatrixComponents%dmu2 = ddualPrimalCoupling * &
!            dtheta * (-p_rspaceTimeDiscr%dgammaC)
!            !dtheta * (-p_rspaceTimeDiscr%dtstep) ! probably wrong!
!
!      END IF

    END IF

  END SUBROUTINE  
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE cc_spaceTimeMatVec (rproblem, rspaceTimeMatrix, rx, rd, cx, cy, &
      cfilter, dnorm, bprintRes)

!<description>
  ! This routine performs a matrix-vector multiplication with the
  ! system matrix A defined by rspaceTimeDiscr.
  !    rd  :=  cx A(p_rsolution) rx  +  cy rd
  ! If rspaceTimeDiscr does not specify an evaluation point for th nonlinearity
  ! in A(.), the routine calculates
  !    rd  :=  cx A(rx) rx  +  cy rd
  ! rd is overwritten by the result.
!</description>

!<input>
  ! A problem structure that provides information about matrices on all
  ! levels as well as temporary vectors.
  TYPE(t_problem), INTENT(INOUT), TARGET :: rproblem

  ! A t_ccoptSpaceTimeDiscretisation structure defining the discretisation of the
  ! coupled space-time matrix.
  TYPE(t_ccoptSpaceTimeMatrix), INTENT(IN) :: rspaceTimeMatrix

  ! A space-time vector defining the current solution.
  TYPE(t_spacetimeVector), INTENT(IN) :: rx
  
  ! Type of filter to apply to the vectors. A combination of SPTID_FILTER_xxxx
  ! flags that specifies which type of filter is to be applied to the
  ! output vector during the matriux vector multiplication.
  INTEGER(I32), INTENT(IN) :: cfilter
  
  ! If set to TRUE and dnorm is present, too, the residuals are printed to the
  ! terminal.
  LOGICAL, INTENT(IN), OPTIONAL :: bprintRes
!</input>

!<inputoutput>
  ! A second space-time vector.
  TYPE(t_spacetimeVector), INTENT(INOUT) :: rd
  
  ! Multiplication factor for rx
  REAL(DP), INTENT(IN) :: cx
  
  ! Multiplication factor for rd
  REAL(DP), INTENT(IN) :: cy
!</inputoutput>

!<output>
  ! OPTIONAL: If specified, returns the $l_2$-norm of rd.
  REAL(DP), INTENT(OUT), OPTIONAL :: dnorm
!<output>

!</subroutine>

    ! local variables
    INTEGER :: ieqTime,ilevel,icp
    TYPE(t_vectorBlock) :: rtempVectorD
    TYPE(t_vectorBlock), DIMENSION(3) :: rtempVector
    TYPE(t_vectorBlock), DIMENSION(3) :: rtempVectorEval
    TYPE(t_blockDiscretisation), POINTER :: p_rdiscr
    REAL(DP) :: dtheta,dnormpart
    TYPE(t_matrixBlock) :: rblockTemp
    TYPE(t_ccmatrixComponents) :: rmatrixComponents
    TYPE(t_ccoptSpaceTimeDiscretisation), POINTER :: p_rspaceTimeDiscretisation
    REAL(DP) :: dtstep
    
    ! DEBUG!!!
    REAL(DP), DIMENSION(:), POINTER :: p_Dx1,p_Dx2,p_Dx3,p_Db
    REAL(DP), DIMENSION(:), POINTER :: p_DxE1,p_DxE2,p_DxE3
    
    ! Pointer to the space time discretisation
    p_rspaceTimeDiscretisation => rspaceTimeMatrix%p_rspaceTimeDiscretisation
    
    ! Level of the discretisation
    ilevel = p_rspaceTimeDiscretisation%ilevel

    ! Theta-scheme identifier.
    ! =1: implicit Euler.
    ! =0.5: Crank Nicolson
    dtheta = rproblem%rtimedependence%dtimeStepTheta
    
    ! Create a temp vector that contains the part of rd which is to be modified.
    p_rdiscr => p_rspaceTimeDiscretisation%p_RlevelInfo%p_rdiscretisation
    CALL lsysbl_createVecBlockByDiscr (p_rdiscr,rtempVectorD,.FALSE.)
    
    ! The vector will be a defect vector. Assign the boundary conditions so
    ! that we can implement them.
    rtempVectorD%p_rdiscreteBC => p_rspaceTimeDiscretisation%p_rlevelInfo%p_rdiscreteBC
    rtempVectorD%p_rdiscreteBCfict => p_rspaceTimeDiscretisation%p_rlevelInfo%p_rdiscreteFBC
    
    ! Create a temp vector for the X-vectors at timestep i-1, i and i+1.
    CALL lsysbl_createVecBlockByDiscr (p_rdiscr,rtempVector(1),.FALSE.)
    CALL lsysbl_createVecBlockByDiscr (p_rdiscr,rtempVector(2),.FALSE.)
    CALL lsysbl_createVecBlockByDiscr (p_rdiscr,rtempVector(3),.FALSE.)
    
    ! Get the parts of the X-vector which are to be modified at first --
    ! subvector 1, 2 and 3. rtempVector(1) contains the 'previous' solution,
    ! rtempVector(2) the 'current' and rtempVector(3) the 'next' one.
    CALL sptivec_getTimestepData(rx, 1+0, rtempVector(2))
      
    ! If necesary, multiply the rtempVectorX. We have to take a -1 into
    ! account as the actual matrix multiplication routine cc_assembleDefect
    ! introduces another -1!
    IF (cx .NE. -1.0_DP) THEN
      CALL lsysbl_scaleVector (rtempVector(2),-cx)
    END IF
      
    ! Now what's the evaluation point where to evaluate the nonlinearity/ies?
    ! If the structure defines an evaluation point, we take that one, so
    ! we need another temp vector that holds the evaluation point.
    ! If not, we take the vector rx. This is archived by creating new
    ! vectors that share their information with rx.
    IF (ASSOCIATED(rspaceTimeMatrix%p_rsolution)) THEN
      CALL lsysbl_createVecBlockByDiscr (p_rdiscr,rtempVectorEval(1),.FALSE.)
      CALL lsysbl_createVecBlockByDiscr (p_rdiscr,rtempVectorEval(2),.FALSE.)
      CALL lsysbl_createVecBlockByDiscr (p_rdiscr,rtempVectorEval(3),.FALSE.)
    ELSE
      CALL lsysbl_duplicateVector (rtempVector(1),rtempVectorEval(1),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      CALL lsysbl_duplicateVector (rtempVector(2),rtempVectorEval(2),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
      CALL lsysbl_duplicateVector (rtempVector(3),rtempVectorEval(3),&
          LSYSSC_DUP_SHARE,LSYSSC_DUP_SHARE)
    END IF

    ! If a nonlinearity is involved, rtempVectorEval(1) contains the 'previous',
    ! rtempVectorEval(2) the 'current' and rtempVectorEval(3) the 'next'
    ! solution (relative to the current time step) where to evaluate
    ! the nonlinearity.
    IF (ASSOCIATED(rspaceTimeMatrix%p_rsolution)) THEN
      CALL sptivec_getTimestepData(rspaceTimeMatrix%p_rsolution, &
          1+0, rtempVectorEval(2))
    END IF
    
    ! DEBUG!!!
    CALL lsysbl_getbase_double (rtempVector(1),p_Dx1)
    CALL lsysbl_getbase_double (rtempVector(2),p_Dx2)
    CALL lsysbl_getbase_double (rtempVector(3),p_Dx3)
    CALL lsysbl_getbase_double (rtempVectorD,p_Db)
    CALL lsysbl_getbase_double (rtempVectorEval(1),p_DxE1)
    CALL lsysbl_getbase_double (rtempVectorEval(2),p_DxE2)
    CALL lsysbl_getbase_double (rtempVectorEval(3),p_DxE3)
    
    ! Basic initialisation of rmatrixComponents with the pointers to the
    ! matrices / discretisation structures on the current level.
    !
    ! The weights in the rmatrixComponents structure are later initialised
    ! according to the actual situation when the matrix is to be used.
    rmatrixComponents%p_rdiscretisation         => &
        p_rspaceTimeDiscretisation%p_rlevelInfo%p_rdiscretisation
    rmatrixComponents%p_rmatrixStokes         => &
        p_rspaceTimeDiscretisation%p_rlevelInfo%rmatrixStokes          
    rmatrixComponents%p_rmatrixB1             => &
        p_rspaceTimeDiscretisation%p_rlevelInfo%rmatrixB1              
    rmatrixComponents%p_rmatrixB2             => &
        p_rspaceTimeDiscretisation%p_rlevelInfo%rmatrixB2              
    rmatrixComponents%p_rmatrixMass           => &
        p_rspaceTimeDiscretisation%p_rlevelInfo%rmatrixMass            
    rmatrixComponents%p_rmatrixIdentityPressure => &
        p_rspaceTimeDiscretisation%p_rlevelInfo%rmatrixIdentityPressure
    rmatrixComponents%iupwind1 = collct_getvalue_int (rproblem%rcollection,'IUPWIND1')
    rmatrixComponents%iupwind2 = collct_getvalue_int (rproblem%rcollection,'IUPWIND2')
    rmatrixComponents%dnu = collct_getvalue_real (rproblem%rcollection,'NU')
    rmatrixComponents%dupsam1 = collct_getvalue_real (rproblem%rcollection,'UPSAM1')
    rmatrixComponents%dupsam2 = collct_getvalue_real (rproblem%rcollection,'UPSAM2')

    ! If dnorm is specified, clear it.
    IF (PRESENT(dnorm)) THEN
      dnorm = 0.0_DP
    END IF

    dtstep = p_rspaceTimeDiscretisation%rtimeDiscr%dtstep
    
    ! Loop through the substeps
    
    DO ieqTime = 0,p_rspaceTimeDiscretisation%NEQtime-1
    
      ! Depending on the time step scheme, initialise the current time.
      SELECT CASE (p_rspaceTimeDiscretisation%rtimeDiscr%ctype)
      CASE (TDISCR_THETA)
        rproblem%rtimedependence%dtime = &
            rproblem%rtimedependence%dtimeInit + ieqTime*dtstep
      CASE (TDISCR_DG0)
        rproblem%rtimedependence%dtime = &
            rproblem%rtimedependence%dtimeInit + (REAL(ieqTime,DP)+0.5_DP)*dtstep
      END SELECT

      ! Get the part of rd which is to be modified.
      IF (cy .NE. 0.0_DP) THEN
        CALL sptivec_getTimestepData(rd, 1+ieqTime, rtempVectorD)
        
        ! If cy <> 1, multiply rtempVectorD by that.
        IF (cy .NE. 1.0_DP) THEN
          CALL lsysbl_scaleVector (rtempVectorD,cy)
        END IF
      ELSE
        CALL lsysbl_clearVector (rtempVectorD)
      END IF
      
      IF (ieqTime .NE. p_rspaceTimeDiscretisation%NEQtime-1) THEN
      
        ! Read the solution of the 'next' timestep and the 'next' evaluation
        ! point.

        CALL sptivec_getTimestepData(rx, 1+ieqTime+1, rtempVector(3))
          
        ! If necesary, multiply the rtempVectorX. We have to take a -1 into
        ! account as the actual matrix multiplication routine cc_assembleDefect
        ! introduces another -1!
        IF (cx .NE. -1.0_DP) THEN
          CALL lsysbl_scaleVector (rtempVector(3),-cx)
        END IF

        IF (ASSOCIATED(rspaceTimeMatrix%p_rsolution)) THEN
          CALL sptivec_getTimestepData(rspaceTimeMatrix%p_rsolution, &
              1+ieqTime+1, rtempVectorEval(3))
        END IF

      END IF

      ! -----
      ! Discretise the boundary conditions at the new point in time -- 
      ! if the boundary conditions are nonconstant in time!
      IF (collct_getvalue_int (rproblem%rcollection,'IBOUNDARY') .NE. 0) THEN
        CALL cc_updateDiscreteBC (rproblem, .FALSE.)
      END IF
      
      ! The first and last substep is a little bit special concerning
      ! the matrix!
      IF (ieqTime .EQ. 0) THEN
        
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
        CALL cc_setupMatrixWeights (rproblem,rspaceTimeMatrix,dtheta,&
          ieqTime,0,rmatrixComponents)
          
        ! Subtract: rd = rd - A11 x1
        CALL cc_assembleDefect (rmatrixComponents,rtempVector(2),rtempVectorD,&
            1.0_DP,rtempVectorEval(1),rtempVectorEval(2),rtempVectorEval(3))

        ! -----
      
        ! Create the matrix
        !   A12  :=  -M + dt*dtheta*[-nu\Laplace u + u \grad u]
        ! and subtract A12 x2 from rd.

        ! Set up the matrix weights of that submatrix.
        CALL cc_setupMatrixWeights (rproblem,rspaceTimeMatrix,dtheta,&
          ieqTime,1,rmatrixComponents)

        ! Subtract: rd = rd - A12 x2
        CALL cc_assembleDefect (rmatrixComponents,rtempVector(3),rtempVectorD,&
            1.0_DP,rtempVectorEval(1),rtempVectorEval(2),rtempVectorEval(3))

        ! Release the block mass matrix.
        CALL lsysbl_releaseMatrix (rblockTemp)

      ELSE IF (ieqTime .LT. p_rspaceTimeDiscretisation%NEQtime-1) THEN

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
        CALL cc_setupMatrixWeights (rproblem,rspaceTimeMatrix,dtheta,&
          ieqTime,-1,rmatrixComponents)
            
        ! Subtract: rd = rd - Aii-1 xi-1.
        ! Note that at this point, the nonlinearity must be evaluated
        ! at xi due to the discretisation scheme!!!
        CALL cc_assembleDefect (rmatrixComponents,rtempVector(1),rtempVectorD,&
            1.0_DP,rtempVectorEval(1),rtempVectorEval(2),rtempVectorEval(3))

        ! Release the block mass matrix.
        CALL lsysbl_releaseMatrix (rblockTemp)

        ! -----      

        ! Now the diagonal matrix.
      
        ! Assemble the nonlinear defect.
      
        ! Set up the matrix weights of that submatrix.
        CALL cc_setupMatrixWeights (rproblem,rspaceTimeMatrix,dtheta,&
          ieqTime,0,rmatrixComponents)

        ! Subtract: rd = rd - Aii xi
        CALL cc_assembleDefect (rmatrixComponents,rtempVector(2),rtempVectorD,&
            1.0_DP,rtempVectorEval(1),rtempVectorEval(2),rtempVectorEval(3))
            
        ! -----
        
        ! Create the matrix
        !   Aii+1 := -M + dt*dtheta*[-nu\Laplace u + u \grad u]
        ! and include that into the global matrix for the dual velocity.

        ! Set up the matrix weights of that submatrix.
        CALL cc_setupMatrixWeights (rproblem,rspaceTimeMatrix,dtheta,&
          ieqTime,1,rmatrixComponents)
          
        ! Subtract: rd = rd - Aii+1 xi+1
        CALL cc_assembleDefect (rmatrixComponents,rtempVector(3),rtempVectorD,&
            1.0_DP,rtempVectorEval(1),rtempVectorEval(2),rtempVectorEval(3))
        
        ! Release the block mass matrix.
        CALL lsysbl_releaseMatrix (rblockTemp)

      ELSE 
        
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
        CALL cc_setupMatrixWeights (rproblem,rspaceTimeMatrix,dtheta,&
          ieqTime,-1,rmatrixComponents)
          
        ! Subtract: rd = rd - Ann-1 xn-1
        ! Note that at this point, the nonlinearity must be evaluated
        ! at xn due to the discretisation scheme!!!
        CALL cc_assembleDefect (rmatrixComponents,rtempVector(1),rtempVectorD,&
            1.0_DP,rtempVectorEval(1),rtempVectorEval(2),rtempVectorEval(3))
     
        ! -----
        
        ! Now the diagonal matrix.
      
        ! Assemble the nonlinear defect.
      
        ! Set up the matrix weights of that submatrix.
        CALL cc_setupMatrixWeights (rproblem,rspaceTimeMatrix,dtheta,&
          ieqTime,0,rmatrixComponents)

        ! Subtract: rd = rd - Ann xn
        CALL cc_assembleDefect (rmatrixComponents,rtempVector(2),rtempVectorD,&
            1.0_DP,rtempVectorEval(1),rtempVectorEval(2),rtempVectorEval(3))
      
      END IF
      
      IF (IAND(cfilter,SPTID_FILTER_BCDEF) .NE. 0) THEN
        ! Implement the boundary conditions into the defect.
        CALL vecfil_discreteBCdef (rtempVectorD)
        CALL vecfil_discreteFBCdef (rtempVectorD)      
      END IF

      IF ((ieqTime .EQ. 0) .AND. (IAND(cfilter,SPTID_FILTER_ICDEF) .NE. 0)) THEN
        ! Implement the initial conditions into the defect.
        CALL tbc_implementInitCondDefSingle (p_rspaceTimeDiscretisation, rtempVectorD)
      END IF

      IF ((ieqTime .EQ. p_rspaceTimeDiscretisation%NEQtime-1) .AND. &
          (IAND(cfilter,SPTID_FILTER_TCDEF) .NE. 0)) THEN
        ! Implement the initial conditions into the defect.
        CALL tbc_implementTermCondDefSingle (p_rspaceTimeDiscretisation, rtempVectorD)
      END IF
      
      ! Save the defect vector back to rd.
      CALL sptivec_setTimestepData(rd, 1+ieqTime, rtempVectorD)
      
      ! If dnorm is specified, calculate the norm of the sub-defect vector and
      ! add it to dnorm.
      IF (PRESENT(dnorm)) THEN
        dnormpart = lsysbl_vectorNorm(rtempVectorD,LINALG_NORML2)**2
        dnorm = dnorm + dnormpart
        
        IF (PRESENT(bprintRes)) THEN
          IF (bprintRes) THEN
            CALL output_line ('||D_'//TRIM(sys_siL(1+ieqTime,10))//'|| = '//&
                TRIM(sys_sdEL(SQRT(dnormpart),10)) )
            DO icp=1,rtempVectorD%nblocks
              CALL output_line ('  ||D_'//&
                  TRIM(sys_siL(1+ieqTime,10))//'^'//TRIM(sys_siL(icp,2))&
                  //'|| = '//&
                  TRIM(sys_sdEL(&
                      lsyssc_vectorNorm(rtempVectorD%RvectorBlock(icp),LINALG_NORML2),10)) )
            END DO
          END IF
        END IF
      END IF
      
      ! Cycle the solution vectors and the evaluation points: 1 <- 2 <- 3.
      CALL lsysbl_copyVector (rtempVector(2),rtempVector(1))
      CALL lsysbl_copyVector (rtempVector(3),rtempVector(2))

      IF (ASSOCIATED(rspaceTimeMatrix%p_rsolution)) THEN
        CALL lsysbl_copyVector (rtempVectorEval(2),rtempVectorEval(1))
        CALL lsysbl_copyVector (rtempVectorEval(3),rtempVectorEval(2))
      END IF
      
      ! Now, the 3rd subvector is again free to take the next vector.
      
    END DO
    
    ! If dnorm is specified, normalise it.
    ! It was calculated from rspaceTimeDiscr%niterations+1 subvectors.
    IF (PRESENT(dnorm)) THEN
      dnorm = SQRT(dnorm / REAL(p_rspaceTimeDiscretisation%NEQtime,DP))
    END IF
    
    ! Release the temp vectors.
    CALL lsysbl_releaseVector (rtempVectorEval(3))
    CALL lsysbl_releaseVector (rtempVectorEval(2))
    CALL lsysbl_releaseVector (rtempVectorEval(1))
    CALL lsysbl_releaseVector (rtempVector(3))
    CALL lsysbl_releaseVector (rtempVector(2))
    CALL lsysbl_releaseVector (rtempVector(1))
    CALL lsysbl_releaseVector (rtempVectorD)
    
  END SUBROUTINE 
   
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
!    TYPE(t_ccmatrixComponents) :: rmatrixComponents
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
!          CALL cc_updateDiscreteBC (rproblem, .FALSE.)
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
  
END MODULE
