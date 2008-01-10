!##############################################################################
!# ****************************************************************************
!# <name> cc2dmediumm2spacetimediscretisation </name>
!# ****************************************************************************
!#
!# <purpose>
!# This module contains routines for the space time discretisation.
!# The structure t_ccoptSpaceTimeDiscretisation defines the shape of the
!# space-time discretisation and describes how to do a matrix-vector
!# multiplication with the global matrix.
!# The routines in this module form the basic set of discretisation routines:
!#
!# 1.) c2d2_initParamsSupersystem
!#     -> Initialises a space-time discretisation according to the parameters
!#        in a DAT file.
!#
!# 2.) c2d2_doneParamsSupersystem
!#     -> Cleans up a space-time discretisation structure allocated with 
!#        c2d2_initParamsSupersystem.
!#
!# 3.) c2d2_spaceTimeMatVec
!#     -> Matrix-Vector multiplication of a space-time vector with the
!#        global matrix.
!#
!# Auxiliary routines:
!#
!# 1.) c2d2_setupMatrixWeights
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
!#                   N = N(y_i) = -nu*Laplace(.) + y_i*grad(.) - (.)*grad(y_i)
!#                   R = R(l_i) = -M                                   
!#
!#                   Newton iteration:
!#                   A = A(y_i) = -nu*Laplace(.) + y_i*grad(.) + grad(y_i)*(.)
!#                   N = N(y_i) = -nu*Laplace(.) + y_i*grad(.) - grad(y_i)*(.)
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
!#                   N(y_i) = -nu*Laplace(.) + y_i*grad(.) - (.)*grad(y_i)
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

MODULE cc2dmediumm2spacetimediscret

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
  USE spacetimevectors
  USE dofmapping
  
  USE matrixio
    
  IMPLICIT NONE

!<types>

!<typeblock>

  ! Defines the basic shape of the supersystem which realises the coupling
  ! between all timesteps.
  TYPE t_ccoptSpaceTimeDiscretisation
  
    ! Spatial refinement level of this matrix.
    INTEGER :: ilevel = 0
  
    ! Number of time steps
    !INTEGER :: niterations         = 0

    ! Absolute start time of the simulation
    !REAL(DP) :: dtimeInit          = 0.0_DP     
    
    ! Maximum time of the simulation
    !REAL(DP) :: dtimeMax           = 1.0_DP
    
    ! Time step length of the time discretisation
    !REAL(DP) :: dtstep

    ! Defines the discretisation in time
    TYPE(t_timeDiscretisation) :: rtimeDiscr

    ! Time-stepping scheme;
    ! 0=one step FD scheme (Euler, CN)
    ! 2=dG(0)
    !INTEGER :: ctimeStepScheme     = 0
    
    ! Parameter for one step scheme (THETA) if itimeStepScheme=0;
    ! =0:Forward Euler(instable), =1: Backward Euler, =0.5: Crank-Nicolson
    !REAL(DP) :: dtimeStepTheta     = 1.0_DP
    
    ! Regularisation parameter for the control $\alpha$. Must be <> 0.
    ! A value of 0.0 disables the terminal condition.
    REAL(DP) :: dalphaC = 1.0_DP
    
    ! Regularisation parameter for the terminal condition 
    ! $\gamma/2*||y(T)-z(T)||$.
    ! A value of 0.0 disables the terminal condition.
    REAL(DP) :: dgammaC = 0.0_DP

    ! Problem-related structure that provides the templates for
    ! matrices/vectors on the spatial level of the matrix.
    TYPE(t_problem_lvl), POINTER :: p_rlevelInfo

  END TYPE

!</typeblock>

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

CONTAINS

  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE c2d2_initParamsSupersystem (rproblem,ilevelTime,ilevelSpace,&
      rspaceTimeDiscr, rx, rb, rd)
  
!<description>
  ! Initialises the time stepping scheme according to the parameters in the
  ! DAT file and the problem structure.
!</description>

!<input>
  ! The problem structure describing the whole problem.
  TYPE(t_problem), INTENT(IN), TARGET :: rproblem
  
  ! 'Refinement level in time'. =1: Initialise as described in
  ! rproblem. >1: Refine (irefineInTime-1) times regularly in time.
  ! (-> #timesteps * 2**(irefineInTime-1) )
  INTEGER, INTENT(IN) :: ilevelTime

  ! 'Refinement level in space'. =1: Calculate on the coarse mesh
  ! >0: calculate on level ilevelSpace. Must be <= rproblem%NLMAX,
  ! >= rproblem%NLMIN.
  INTEGER, INTENT(IN) :: ilevelSpace
!</input>

!<inputoutput>
  ! Supersystem-structure to be initialised.
  TYPE(t_ccoptSpaceTimeDiscretisation), INTENT(OUT) :: rspaceTimeDiscr

  ! A space-time vector that is initialised for the current solution.
  TYPE(t_spacetimeVector), INTENT(INOUT), OPTIONAL :: rx

  ! A space-time vector that is initialised for the current RHS.
  TYPE(t_spacetimeVector), INTENT(INOUT), OPTIONAL :: rb

  ! A space-time vector that is initialised for the current defect.
  TYPE(t_spacetimeVector), INTENT(INOUT), OPTIONAL :: rd
!</inputoutput>

!</subroutine>

    ! local variables
    INTEGER :: niterations

    ! Copy most relevant data from the problem structure.
    rspaceTimeDiscr%ilevel          = ilevelSpace
    rspaceTimeDiscr%p_rlevelInfo    => rproblem%RlevelInfo(ilevelSpace)
    niterations = rproblem%rtimedependence%niterations * 2**MAX(0,ilevelTime-1)
    
    ! Initialise the time discretisation
    SELECT CASE (rproblem%rtimedependence%ctimeStepScheme)
    CASE (0) 
      CALL tdiscr_initTheta(&
          rproblem%rtimedependence%dtimeInit,&
          rproblem%rtimedependence%dtimeMax,&
          niterations,rproblem%rtimedependence%dtimeStepTheta,&
          rspaceTimeDiscr%rtimeDiscr)
    CASE (2)
      CALL tdiscr_initdG0(rproblem%rtimedependence%dtimeInit,&
          rproblem%rtimedependence%dtimeMax,&
          niterations, rspaceTimeDiscr%rtimeDiscr)
    CASE DEFAULT
      PRINT *,'c2d2_initParamsSupersystem: Unsupported time discretisation.'
      CALL sys_halt()
    END SELECT
    
    CALL parlst_getvalue_double (rproblem%rparamList,'OPTIMALCONTROL',&
                                'dalphaC',rspaceTimeDiscr%dalphaC,1.0_DP)
    CALL parlst_getvalue_double (rproblem%rparamList,'OPTIMALCONTROL',&
                                'dgammaC',rspaceTimeDiscr%dgammaC,0.0_DP)

    ! Initialise the global solution- and defect- vector.
    IF (PRESENT(rx)) THEN
      CALL sptivec_initVector (rx,&
          rspaceTimeDiscr%rtimeDiscr,rspaceTimeDiscr%p_rlevelInfo%p_rdiscretisation)
    END IF
    IF (PRESENT(rb)) THEN
      CALL sptivec_initVector (rb,&
          rspaceTimeDiscr%rtimeDiscr,rspaceTimeDiscr%p_rlevelInfo%p_rdiscretisation)
    END IF
    IF (PRESENT(rd)) THEN
      CALL sptivec_initVector (rd,&
          rspaceTimeDiscr%rtimeDiscr,rspaceTimeDiscr%p_rlevelInfo%p_rdiscretisation)
    END IF

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE c2d2_doneParamsSupersystem (rspaceTimeDiscr,rx,rb,rd)
  
!<description>
  ! Cleans up a given supersystem structure.
!</description>

!<inputoutput>
  ! Supersystem-structure to be cleaned up.
  TYPE(t_ccoptSpaceTimeDiscretisation), INTENT(OUT) :: rspaceTimeDiscr

  ! A space-time vector that is initialised for the current solution.
  TYPE(t_spacetimeVector), INTENT(INOUT), OPTIONAL :: rx

  ! A space-time vector that is initialised for the current RHS.
  TYPE(t_spacetimeVector), INTENT(INOUT), OPTIONAL :: rb

  ! A space-time vector that is initialised for the current defect.
  TYPE(t_spacetimeVector), INTENT(INOUT), OPTIONAL :: rd
!</inputoutput>

!</subroutine>

    ! Release memory.
    IF (PRESENT(rb)) CALL sptivec_releaseVector (rb)
    IF (PRESENT(rd)) CALL sptivec_releaseVector (rd)
    IF (PRESENT(rx)) CALL sptivec_releaseVector (rx)

  END SUBROUTINE

  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE c2d2_setupMatrixWeights (rproblem,rspaceTimeMatrix,dtheta,&
      isubstep,irelpos,rmatrixComponents,ivecIndexNonlin)

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
  ! A t_ccmatrixComponents structure that defines tghe shape of the core
  ! equation. The weights that specify the submatrices of a small 6x6 
  ! block matrix system are initialised depending on the position
  ! specified by isubstep and nsubsteps.
  TYPE(t_ccmatrixComponents), INTENT(INOUT) :: rmatrixComponents
  
  ! Returns the (1-based) number of the subvector in rspaceTimeMatrix%p_rsolution
  ! which should be specified as evaluation point for the nonlinearity.
  ! If no nonlinearity has to be evaluated, a negative value is returned.
  INTEGER, INTENT(OUT) :: ivecIndexNonlin
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
    
    p_rspaceTimeDiscr => rspaceTimeMatrix%p_rspaceTimeDiscretisation
    
    dequationType = 1.0_DP
    IF (rproblem%roptcontrol%ispaceTimeFormulation .NE. 0) &
      dequationType = -1.0_DP

    ! What's the matrix type we should set up? If we have to set up a
    ! Newton matrix, we put dnewton to 1.0 so the Newton part in the
    ! primal velocity is assembled.
    dnewton = 0.0_DP
    IF (rspaceTimeMatrix%cmatrixType .EQ. 1) THEN
      IF (rproblem%iequation .EQ. 0) THEN
        ! Newton is only to be assembled in Navier-Stokes!
        dnewton = 1.0_DP
      END IF
      dsmooth = 1.0_DP
    ELSE 
      ! Activate smoothing in the last timestep.
      dsmooth = 1.0_DP
    END IF

    dtstep = p_rspaceTimeDiscr%rtimeDiscr%dtstep

    ! The first and last substep is a little bit special concerning
    ! the matrix!
    IF (isubstep .EQ. 0) THEN
      
      ! We are in the first substep
    
      IF (irelpos .EQ. 0) THEN
      
        ! The diagonal matrix.
        ivecIndexNonlin = 1 + isubstep
        
        rmatrixComponents%diota1 = 1.0_DP
        rmatrixComponents%diota2 = 0.0_DP

        rmatrixComponents%dkappa1 = 1.0_DP
        rmatrixComponents%dkappa2 = 0.0_DP
        
        rmatrixComponents%dalpha1 = 0.0_DP
        rmatrixComponents%dalpha2 = dtimeCoupling * 1.0_DP/dtstep
        
        rmatrixComponents%dtheta1 = 0.0_DP
        rmatrixComponents%dtheta2 = dtheta
        
        rmatrixComponents%dgamma1 = 0.0_DP
        rmatrixComponents%dgamma2 = &
            - dtheta * REAL(1-rproblem%iequation,DP)
        
        rmatrixComponents%dnewton1 = 0.0_DP
        rmatrixComponents%dnewton2 = &
              dtheta * REAL(1-rproblem%iequation,DP)

        rmatrixComponents%deta1 = 0.0_DP
        rmatrixComponents%deta2 = 1.0_DP
        
        rmatrixComponents%dtau1 = 0.0_DP
        rmatrixComponents%dtau2 = 1.0_DP
        
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
          rmatrixComponents%dmu2 = ddualPrimalCoupling * &
              dtheta * (-1.0_DP)
          rmatrixComponents%dr21 = 0.0_DP
          rmatrixComponents%dr22 = 0.0_DP
        ELSE
          rmatrixComponents%dmu2 = ddualPrimalCoupling * &
              dtheta * (-1.0_DP)
          rmatrixComponents%dr21 = ddualPrimalCoupling * &
              dtheta * ( 1.0_DP)
          rmatrixComponents%dr22 = ddualPrimalCoupling * &
              dtheta * (-1.0_DP)
        END IF
        rmatrixComponents%dmu1 = 0.0_DP
        !rmatrixComponents%dmu2 = 0.0_DP
                    
      ELSE IF (irelpos .EQ. 1) THEN
      
        ! Offdiagonal matrix on the right of the diagonal.
        ivecIndexNonlin = 1 + isubstep+1

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
        
        rmatrixComponents%dgamma1 = 0.0_DP
        rmatrixComponents%dgamma2 = &
            - (1.0_DP-dtheta) * REAL(1-rproblem%iequation,DP)
        
        rmatrixComponents%dnewton1 = 0.0_DP
        rmatrixComponents%dnewton2 = &
              (1.0_DP-dtheta) * REAL(1-rproblem%iequation,DP)

        rmatrixComponents%deta1 = 0.0_DP
        rmatrixComponents%deta2 = 0.0_DP
        
        rmatrixComponents%dtau1 = 0.0_DP
        rmatrixComponents%dtau2 = 0.0_DP
        
        rmatrixComponents%dmu1 = 0.0_DP
        rmatrixComponents%dmu2 = ddualPrimalCoupling * &
            (-dequationType) * (1.0_DP-dtheta)
            
        rmatrixComponents%dr21 = 0.0_DP
        rmatrixComponents%dr22 = 0.0_DP
            
      END IF
    
    ELSE IF (isubstep .LT. p_rspaceTimeDiscr%rtimeDiscr%nintervals) THEN
      
      ! We are sonewhere in the middle of the matrix. There is a substep
      ! isubstep+1 and a substep isubstep-1!
      
      ! -----
      
      IF (irelpos .EQ. -1) THEN
      
        ! Matrix on the left of the diagonal.
        !
        ! Note that at this point, the nonlinearity must be evaluated
        ! at xi due to the discretisation scheme!!!
        ivecIndexNonlin = 1 + isubstep

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
        
        rmatrixComponents%dgamma1 = &
            (1.0_DP-dtheta) * REAL(1-rproblem%iequation,DP)
        rmatrixComponents%dgamma2 = 0.0_DP
        
        rmatrixComponents%dnewton1 = 0.0_DP
        rmatrixComponents%dnewton2 = 0.0_DP

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
        ivecIndexNonlin = 1 + isubstep

        rmatrixComponents%diota1 = 0.0_DP
        rmatrixComponents%diota2 = 0.0_DP

        rmatrixComponents%dkappa1 = 0.0_DP
        rmatrixComponents%dkappa2 = 0.0_DP
        
        rmatrixComponents%dalpha1 = dtimeCoupling * 1.0_DP/dtstep
        rmatrixComponents%dalpha2 = dtimeCoupling * 1.0_DP/dtstep
        
        rmatrixComponents%dtheta1 = dtheta
        rmatrixComponents%dtheta2 = dtheta
        
        rmatrixComponents%dgamma1 = &
            dtheta * REAL(1-rproblem%iequation,DP)
        rmatrixComponents%dgamma2 = &
            - dtheta * REAL(1-rproblem%iequation,DP)
        
        rmatrixComponents%dnewton1 = dtheta * dnewton
        rmatrixComponents%dnewton2 = &
              dtheta * REAL(1-rproblem%iequation,DP)

        rmatrixComponents%deta1 = 1.0_DP
        rmatrixComponents%deta2 = 1.0_DP
        
        rmatrixComponents%dtau1 = 1.0_DP
        rmatrixComponents%dtau2 = 1.0_DP
        
        rmatrixComponents%dmu1 = dprimalDualCoupling * &
            dequationType * dtheta * 1.0_DP / p_rspaceTimeDiscr%dalphaC
        IF (dnewton .EQ. 0.0_DP) THEN
          rmatrixComponents%dmu2 = ddualPrimalCoupling * &
              (-dequationType) * dtheta 
          rmatrixComponents%dr21 = 0.0_DP
          rmatrixComponents%dr22 = 0.0_DP
        ELSE
          rmatrixComponents%dmu2 = ddualPrimalCoupling * &
              (-dequationType) * dtheta 
          rmatrixComponents%dr21 = ddualPrimalCoupling * &
              ( dequationType) * dtheta 
          rmatrixComponents%dr22 = ddualPrimalCoupling * &
              (-dequationType) * dtheta 
        END IF

      ELSE IF (irelpos .EQ. 1) THEN
            
        ! Matrix on the right of the diagonal.
        ivecIndexNonlin = 1 + isubstep+1

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
        
        rmatrixComponents%dgamma1 = 0.0_DP
        rmatrixComponents%dgamma2 = &
            - (1.0_DP-dtheta) * REAL(1-rproblem%iequation,DP)
        
        rmatrixComponents%dnewton1 = 0.0_DP
        rmatrixComponents%dnewton2 = &
            (1.0_DP-dtheta) * REAL(1-rproblem%iequation,DP)

        rmatrixComponents%deta1 = 0.0_DP
        rmatrixComponents%deta2 = 0.0_DP
        
        rmatrixComponents%dtau1 = 0.0_DP
        rmatrixComponents%dtau2 = 0.0_DP
        
        rmatrixComponents%dmu1 = 0.0_DP
        rmatrixComponents%dmu2 = ddualPrimalCoupling * &
            (-dequationType) * (1.0_DP-dtheta) 
            
        rmatrixComponents%dr21 = 0.0_DP
        rmatrixComponents%dr22 = 0.0_DP
            
      END IF
    
    ELSE
    
      ! We are in the last substep
      
!      ! Although this would be the correct implementation of the matrix weights
!      ! for the terminal condition at the first glance, it would be wrong to use 
!      ! that. The last timestep has to be processed like
!      ! the others, which infers some 'smoothing' in the dual solution
!      ! at the end of the time cylinder! Without that smoothing that, some 
!      ! time discretisation schemes like Crank-Nicolson would get instable
!      ! if a terminal condition <> 0 is prescribed!
!      
!      IF (irelpos .EQ. -1) THEN
!      
!        ! Matrix on the left of the diagonal
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
!        rmatrixComponents%diota2 = (1.0_DP-ddualPrimalCoupling)
!
!        rmatrixComponents%dkappa1 = 0.0_DP
!        rmatrixComponents%dkappa2 = 1.0_DP
!        
!        rmatrixComponents%dalpha1 = 1.0_DP
!        rmatrixComponents%dalpha2 = 1.0_DP
!        
!        rmatrixComponents%dtheta1 = dtheta * p_rspaceTimeDiscr%dtstep
!        rmatrixComponents%dtheta2 = 0.0_DP
!        
!        rmatrixComponents%dgamma1 = &
!            dtheta * p_rspaceTimeDiscr%dtstep * REAL(1-rproblem%iequation,DP)
!        rmatrixComponents%dgamma2 = 0.0_DP
!!           - dtheta * p_rspaceTimeDiscr%dtstep * REAL(1-rproblem%iequation,DP)
!        
!        rmatrixComponents%dnewton1 = 0.0_DP
!        rmatrixComponents%dnewton2 = 0.0_DP
!!             dtheta * p_rspaceTimeDiscr%dtstep * REAL(1-rproblem%iequation,DP)
!
!        rmatrixComponents%deta1 = p_rspaceTimeDiscr%dtstep
!        rmatrixComponents%deta2 = 0.0_DP
!        
!        rmatrixComponents%dtau1 = 1.0_DP
!        rmatrixComponents%dtau2 = 0.0_DP
!        
!        rmatrixComponents%dmu1 = dprimalDualCoupling * &
!            dtheta * p_rspaceTimeDiscr%dtstep / p_rspaceTimeDiscr%dalphaC
!        rmatrixComponents%dmu2 = ddualPrimalCoupling * &
!            (-p_rspaceTimeDiscr%dgammaC)
!
!      END IF
        
      IF (irelpos .EQ. -1) THEN
      
        ! Matrix on the left of the diagonal.
        !
        ! Note that at this point, the nonlinearity must be evaluated
        ! at xn due to the discretisation scheme!!!
        ivecIndexNonlin = 1 + isubstep

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
        
        rmatrixComponents%dgamma1 = &
            (1.0_DP-dtheta) * REAL(1-rproblem%iequation,DP)
        rmatrixComponents%dgamma2 = 0.0_DP
        
        rmatrixComponents%dnewton1 = 0.0_DP
        rmatrixComponents%dnewton2 = 0.0_DP

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
        ivecIndexNonlin = 1 + isubstep

        rmatrixComponents%diota1 = 0.0_DP
        rmatrixComponents%diota2 = (1.0_DP-dsmooth) !0.0_DP

        rmatrixComponents%dkappa1 = 0.0_DP
        rmatrixComponents%dkappa2 = (1.0_DP-dsmooth) !0.0_DP
        
        rmatrixComponents%dalpha1 = dtimeCoupling * 1.0_DP/dtstep
        rmatrixComponents%dalpha2 = (1.0_DP-dsmooth) + &
            dsmooth*dtimeCoupling * 1.0_DP/dtstep
        
        rmatrixComponents%dtheta1 = dtheta
        rmatrixComponents%dtheta2 = dsmooth * dtheta
        
        rmatrixComponents%dgamma1 = &
            dtheta * REAL(1-rproblem%iequation,DP)
        rmatrixComponents%dgamma2 = &
            - dsmooth * dtheta * REAL(1-rproblem%iequation,DP)
        
        rmatrixComponents%dnewton1 = dtheta * dnewton
        rmatrixComponents%dnewton2 = &
              dsmooth * dtheta * REAL(1-rproblem%iequation,DP)

        rmatrixComponents%deta1 = 1.0_DP
        rmatrixComponents%deta2 = dsmooth
        
        rmatrixComponents%dtau1 = 1.0_DP
        rmatrixComponents%dtau2 = dsmooth
        
        rmatrixComponents%dmu1 = dprimalDualCoupling * &
            dequationType * dtheta * 1.0_DP / p_rspaceTimeDiscr%dalphaC
            
        ! Weight the mass matrix by GAMMA instead of delta(T).
        ! That's the only difference to the implementation above!
        IF (dnewton .EQ. 0.0_DP) THEN
          rmatrixComponents%dmu2 = ddualPrimalCoupling * &
              (-dequationType) * dtheta * p_rspaceTimeDiscr%dgammaC
              !dtheta * p_rspaceTimeDiscr%dtstep
          rmatrixComponents%dr21 = 0.0_DP
          rmatrixComponents%dr22 = 0.0_DP
        ELSE
          rmatrixComponents%dmu2 = ddualPrimalCoupling * &
              (-dequationType) * dtheta * p_rspaceTimeDiscr%dgammaC
          rmatrixComponents%dr21 = ddualPrimalCoupling * &
              ( dequationType) * dtheta 
          rmatrixComponents%dr22 = ddualPrimalCoupling * &
              (-dequationType) * dtheta 
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

    ! If there is no nonlinearity involved, return 0 as subvector index 
    IF ((rmatrixComponents%dgamma1 .EQ. 0.0_DP) .AND. &
        (rmatrixComponents%dgamma2 .EQ. 0.0_DP) .AND. &
        (rmatrixComponents%dnewton1 .EQ. 0.0_DP) .AND. &
        (rmatrixComponents%dnewton2 .EQ. 0.0_DP))  THEN
      ivecIndexNonlin = -ivecIndexNonlin
    END IF

  END SUBROUTINE  
  
  ! ***************************************************************************
  
!<subroutine>

  SUBROUTINE c2d2_spaceTimeMatVec (rproblem, rspaceTimeMatrix, rx, rd, cx, cy, dnorm,&
      bprintRes)

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
    INTEGER :: isubstep,ilevel,icp,iidxNonlin,irelNonl
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
    ! subvector 1, 2 and 3.
    CALL sptivec_getTimestepData(rx, 0, rtempVector(1))
    IF (p_rspaceTimeDiscretisation%rtimeDiscr%nintervals .GT. 0) &
      CALL sptivec_getTimestepData(rx, 1, rtempVector(2)) 
    IF (p_rspaceTimeDiscretisation%rtimeDiscr%nintervals .GT. 1) THEN
      CALL sptivec_getTimestepData(rx, 2, rtempVector(3))
    ELSE
      CALL lsysbl_copyVector (rtempVector(2),rtempVector(3))
    END IF
      
    ! If necesary, multiply the rtempVectorX. We have to take a -1 into
    ! account as the actual matrix multiplication routine c2d2_assembleDefect
    ! introduces another -1!
    IF (cx .NE. -1.0_DP) THEN
      CALL lsysbl_scaleVector (rtempVector(1),-cx)
      CALL lsysbl_scaleVector (rtempVector(2),-cx)
      CALL lsysbl_scaleVector (rtempVector(3),-cx)
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

    IF (ASSOCIATED(rspaceTimeMatrix%p_rsolution)) THEN
      CALL sptivec_getTimestepData(rspaceTimeMatrix%p_rsolution, 0, rtempVectorEval(1))
      IF (p_rspaceTimeDiscretisation%rtimeDiscr%nintervals .GT. 0) &
        CALL sptivec_getTimestepData(rspaceTimeMatrix%p_rsolution, 1, rtempVectorEval(2))
      IF (p_rspaceTimeDiscretisation%rtimeDiscr%nintervals .GT. 1) THEN
        CALL sptivec_getTimestepData(rspaceTimeMatrix%p_rsolution, 2, rtempVectorEval(3))
      ELSE
        CALL lsysbl_copyVector (rtempVectorEval(2),rtempVectorEval(3))
      END IF
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
    
    DO isubstep = 0,p_rspaceTimeDiscretisation%rtimeDiscr%nintervals
    
      ! Current point in time
      rproblem%rtimedependence%dtime = &
          rproblem%rtimedependence%dtimeInit + isubstep*dtstep

      ! Get the part of rd which is to be modified.
      IF (cy .NE. 0.0_DP) THEN
        CALL sptivec_getTimestepData(rd, isubstep, rtempVectorD)
        
        ! If cy <> 1, multiply rtempVectorD by that.
        IF (cy .NE. 1.0_DP) THEN
          CALL lsysbl_scaleVector (rtempVectorD,cy)
        END IF
      ELSE
        CALL lsysbl_clearVector (rtempVectorD)
      END IF

      ! -----
      ! Discretise the boundary conditions at the new point in time -- 
      ! if the boundary conditions are nonconstant in time!
      IF (collct_getvalue_int (rproblem%rcollection,'IBOUNDARY') .NE. 0) THEN
        CALL c2d2_updateDiscreteBC (rproblem, .FALSE.)
      END IF
      
      ! The first and last substep is a little bit special concerning
      ! the matrix!
      IF (isubstep .EQ. 0) THEN
        
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
        CALL c2d2_setupMatrixWeights (rproblem,rspaceTimeMatrix,dtheta,&
          isubstep,0,rmatrixComponents,iidxNonlin)
          
        ! Subtract: rd = rd - A11 x1
        irelnonl = ABS(iidxNonlin)-isubstep
        CALL c2d2_assembleDefect (rmatrixComponents,rtempVector(1),rtempVectorD,&
            1.0_DP,rtempVectorEval(irelnonl))

        ! -----
      
        ! Create the matrix
        !   A12  :=  -M + dt*dtheta*[-nu\Laplace u + u \grad u]
        ! and subtract A12 x2 from rd.

        ! Set up the matrix weights of that submatrix.
        CALL c2d2_setupMatrixWeights (rproblem,rspaceTimeMatrix,dtheta,&
          isubstep,1,rmatrixComponents,iidxNonlin)

        ! Subtract: rd = rd - A12 x2
        irelnonl = ABS(iidxNonlin)-isubstep
        CALL c2d2_assembleDefect (rmatrixComponents,rtempVector(2),rtempVectorD,&
            1.0_DP,rtempVectorEval(irelnonl))

        ! Release the block mass matrix.
        CALL lsysbl_releaseMatrix (rblockTemp)

      ELSE IF (isubstep .LT. p_rspaceTimeDiscretisation%rtimeDiscr%nintervals) THEN

        ! We are sonewhere in the middle of the matrix. There is a substep
        ! isubstep+1 and a substep isubstep-1!  Here, we have to handle the following
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
        CALL c2d2_setupMatrixWeights (rproblem,rspaceTimeMatrix,dtheta,&
          isubstep,-1,rmatrixComponents,iidxNonlin)
            
        ! Subtract: rd = rd - Aii-1 xi-1.
        ! Note that at this point, the nonlinearity must be evaluated
        ! at xi due to the discretisation scheme!!!
        irelnonl = ABS(iidxNonlin)-isubstep + 1
        CALL c2d2_assembleDefect (rmatrixComponents,rtempVector(1),rtempVectorD,&
            1.0_DP,rtempVectorEval(irelnonl))

        ! Release the block mass matrix.
        CALL lsysbl_releaseMatrix (rblockTemp)

        ! -----      

        ! Now the diagonal matrix.
      
        ! Assemble the nonlinear defect.
      
        ! Set up the matrix weights of that submatrix.
        CALL c2d2_setupMatrixWeights (rproblem,rspaceTimeMatrix,dtheta,&
          isubstep,0,rmatrixComponents,iidxNonlin)

        ! Subtract: rd = rd - Aii xi
        irelnonl = ABS(iidxNonlin)-isubstep + 1
        CALL c2d2_assembleDefect (rmatrixComponents,rtempVector(2),rtempVectorD,&
            1.0_DP,rtempVectorEval(irelnonl))
            
        ! -----
        
        ! Create the matrix
        !   Aii+1 := -M + dt*dtheta*[-nu\Laplace u + u \grad u]
        ! and include that into the global matrix for the dual velocity.

        ! Set up the matrix weights of that submatrix.
        CALL c2d2_setupMatrixWeights (rproblem,rspaceTimeMatrix,dtheta,&
          isubstep,1,rmatrixComponents,iidxNonlin)
          
        ! Subtract: rd = rd - Aii+1 xi+1
        irelnonl = ABS(iidxNonlin)-isubstep + 1
        CALL c2d2_assembleDefect (rmatrixComponents,rtempVector(3),rtempVectorD,&
            1.0_DP,rtempVectorEval(irelnonl))
        
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
        CALL c2d2_setupMatrixWeights (rproblem,rspaceTimeMatrix,dtheta,&
          isubstep,-1,rmatrixComponents,iidxNonlin)
          
        ! Subtract: rd = rd - Ann-1 xn-1
        ! Note that at this point, the nonlinearity must be evaluated
        ! at xn due to the discretisation scheme!!!
        irelnonl = ABS(iidxNonlin)-isubstep + 2
        CALL c2d2_assembleDefect (rmatrixComponents,rtempVector(2),rtempVectorD,&
            1.0_DP,rtempVectorEval(irelnonl))
     
        ! -----
        
        ! Now the diagonal matrix.
      
        ! Assemble the nonlinear defect.
      
        ! Set up the matrix weights of that submatrix.
        CALL c2d2_setupMatrixWeights (rproblem,rspaceTimeMatrix,dtheta,&
          isubstep,0,rmatrixComponents,iidxNonlin)

        ! Subtract: rd = rd - Ann xn
        irelnonl = ABS(iidxNonlin)-isubstep + 2
        CALL c2d2_assembleDefect (rmatrixComponents,rtempVector(3),rtempVectorD,&
            1.0_DP,rtempVectorEval(irelnonl))
      
      END IF
      
      ! Implement the boundary conditions into the defect.
      CALL vecfil_discreteBCdef (rtempVectorD)
      CALL vecfil_discreteFBCdef (rtempVectorD)      
      
      ! Save the defect vector back to rd.
      CALL sptivec_setTimestepData(rd, isubstep, rtempVectorD)
      
      ! If dnorm is specified, calculate the norm of the sub-defect vector and
      ! add it to dnorm.
      IF (PRESENT(dnorm)) THEN
        dnormpart = lsysbl_vectorNorm(rtempVectorD,LINALG_NORML2)**2
        dnorm = dnorm + dnormpart
        
        IF (PRESENT(bprintRes)) THEN
          IF (bprintRes) THEN
            CALL output_line ('||D_'//TRIM(sys_siL(isubstep,10))//'|| = '//&
                TRIM(sys_sdEL(&
                    SQRT(lsysbl_vectorNorm(rtempVectorD,LINALG_NORML2)&
                    ),10)) )
            DO icp=1,rtempVectorD%nblocks
              CALL output_line ('  ||D_'//TRIM(sys_siL(isubstep,10))//'^'//TRIM(sys_siL(icp,2))&
                  //'|| = '//&
                  TRIM(sys_sdEL(&
                      SQRT(lsyssc_vectorNorm(rtempVectorD%RvectorBlock(icp),LINALG_NORML2)&
                      ),10)) )
            END DO
          END IF
        END IF
      END IF
      
      IF ((isubstep .GT. 0) .AND. &
          (isubstep .LT. p_rspaceTimeDiscretisation%rtimeDiscr%nintervals-1)) THEN
      
        ! Shift the timestep data: x_n+1 -> x_n -> x_n-1
        IF (rtempVector(2)%NEQ .NE. 0) &
          CALL lsysbl_copyVector (rtempVector(2), rtempVector(1))
        
        IF (rtempVector(3)%NEQ .NE. 0) &
          CALL lsysbl_copyVector (rtempVector(3), rtempVector(2))
        
        ! Get the new x_n+1 for the next pass through the loop.
        CALL sptivec_getTimestepData(rx, isubstep+2, rtempVector(3))
        
        ! If necessary, multiply rtempVector(3)
        IF (cx .NE. -1.0_DP) THEN
          CALL lsysbl_scaleVector (rtempVector(3),-cx)
        END IF        
        
        ! A similar shifting has to be done for the evaluation point --
        ! if it's different from rx!
        IF (ASSOCIATED(rspaceTimeMatrix%p_rsolution)) THEN
          CALL lsysbl_copyVector (rtempVectorEval(2), rtempVectorEval(1))
          CALL lsysbl_copyVector (rtempVectorEval(3), rtempVectorEval(2))
          ! Get the new x_n+1 for the next pass through the loop.
          CALL sptivec_getTimestepData(rspaceTimeMatrix%p_rsolution, &
              isubstep+2, rtempVectorEval(3))
        END IF
        
      ELSE IF ((p_rspaceTimeDiscretisation%rtimeDiscr%nintervals .EQ. 1) .AND. &
               (isubstep .EQ. 0)) THEN
      
        ! Only one timestep. Copy vector-1 -> vector-2 -> vector-3 so that the above
        ! assembly routine for the last timestep will work.
        CALL lsysbl_copyVector (rtempVector(2),rtempVector(3))
        CALL lsysbl_copyVector (rtempVectorEval(2),rtempVectorEval(3))

        CALL lsysbl_copyVector (rtempVector(1),rtempVector(2))
        CALL lsysbl_copyVector (rtempVectorEval(1),rtempVectorEval(2))
        
      END IF
    
    END DO
    
    ! If dnorm is specified, normalise it.
    ! It was calculated from rspaceTimeDiscr%niterations+1 subvectors.
    IF (PRESENT(dnorm)) THEN
      dnorm = SQRT(dnorm / REAL(p_rspaceTimeDiscretisation%rtimeDiscr%nintervals+1,DP))
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
!  SUBROUTINE c2d2_assembleSpaceTimeRHS (rproblem, rspaceTimeDiscr, rb, &
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
!        CALL c2d2_generateBasicRHS (rproblem,rvector)
!
!      END IF
!      
!      ! Include BC's?
!      IF (bincludeBC) THEN
!        
!        ! Initialise the collection for the assembly process with callback routines.
!        ! Basically, this stores the simulation time in the collection if the
!        ! simulation is nonstationary.
!        CALL c2d2_initCollectForAssembly (rproblem,rproblem%rcollection)
!
!        ! Discretise the boundary conditions at the new point in time -- 
!        ! if the boundary conditions are nonconstant in time!
!        IF (collct_getvalue_int (rproblem%rcollection,'IBOUNDARY') .NE. 0) THEN
!          CALL c2d2_updateDiscreteBC (rproblem, .FALSE.)
!        END IF
!
!        ! Implement the boundary conditions into the RHS.
!        ! This is done *after* multiplying -z by GAMMA or dtstep, resp.,
!        ! as Dirichlet values mustn't be multiplied with GAMMA!
!        CALL vecfil_discreteBCsol (rvector)
!        CALL vecfil_discreteFBCsol (rvector)      
!      
!        ! Clean up the collection (as we are done with the assembly, that's it.
!        CALL c2d2_doneCollectForAssembly (rproblem,rproblem%rcollection)
!        
!      END IF
!    
!    END SUBROUTINE
!    
!  END SUBROUTINE
  
END MODULE
