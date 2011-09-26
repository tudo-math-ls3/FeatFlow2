************************************************************************
* This file describes a structure block that sets up the
* adaptive time stepping strategy.
* 
* There are integer and double precision variables that configure
* the time stepping. Corresponding to that, there is an integer
* and a double precision structure that contains this configuration.
*
* The structures themselves are realized as integer/double precision
* arrays:
*   IPARAM = array [1..SZADTI] of integer 
*          = TAdaptiveTimeStepIParams
*   DPARAM = array [1..SZADTD] of double
*          = TAdaptiveTimeStepDParams
*
* The whole structure can be initializes with the INIADT subroutine
* in NONSTEAD2.F.
*
* TAdaptiveTimeStepIParams = record 
*   IADTIM : Identifier for adaptive time stepping
*            =0: no adaptive time stepping
*            =1: adaptive time stepping with prediction,
*                no repetition except for when the solver breaks down
*            =2: adaptive time stepping with prediction,
*                repetition of the time step if nonlinear stopping
*                criterion is too large or the solver breaks down
*            =3: adaptive time stepping with prediction,
*                repetition of the time step if nonlinear stopping
*                criterion is too large or the time error is too large
*                or the solver breaks down
*
*   IEXTIM : Calculate a solution with extrapolation in time.
*            =0: No extrapolation in time
*            =1: Use information of predictor and corrector step to
*                extrapolate in time
*            Can only be used if IADTIM <> 0.
*
*   IADIN  : Parameter for error control in the start phase
*            =0: Use EPSADI for controlling the error in the start phase
*                (Standard)
*            =1: Use linear combination of EPSADI and EPSADL for
*                controlling the error in the start phase
*            =2: Use logarithmic combination of EPSADI and EPSADL for
*                controlling the error in the start phase
*
*   IREPIT : Maximum number of repetitions if IADTIM=2,3
*
*   IEPSAD : Type of error control in space/time
*            =1: Calculate RELT=RELU2
*            =2: Calculate RELT=RELUM
*            =3: Calculate RELT=RELP2
*            =4: Calculate RELT=RELPM
*            =5: Calculate RELT=max(RELU2,RELP2)
*            =6: Calculate RELT=max(RELUM,RELPM)
*            =7: Calculate RELT=max(RELU2,RELP2,RELUM,RELPM)
*            =8: Calculate RELT=min(RELU2,RELP2,RELUM,RELPM)
* end
*
* TAdaptiveTimeStepDParams = record
*   DTMIN  : Minimum time step length
*   DTMAX  : Maximum time step length
*   DTFACT : Factor for modifying time step size. 
*            Only if IADTIM<>0.
*            Time step size is reduced by sqrt(DTFACT) upon broke down
*            of nonlinear solver, and by DTFACT upon broke down of
*            linear solver
*   TIMSTP : Initial time step length
*   TIMEIN : Absolute time for start procedure;
*            Time step will not be reduced until TIMEIN is reached.
*            This is essential as Navier Stokes might otherwise
*            stagnate and reduce the time step to 0 at the beginning
*            of the simulation!
*   EPSADI : Parameter for time error limit in the start phase; 
*            only respected if IADTIM=3
*   EPSADL : Parameter for time error limit after the start phase; 
*            only respected if IADTIM=3
*   EPSADU : Upper limit for relation of calculated step sizes to accept
*            a time step; if the relation 
*                     (new step size)/(old step size)
*            is too large, repeat the step;
*            only respected if IADTIM=3
*   EPSNS  : Lower limit for the time derivative to be treated as zero
*            -> solution got stationary
* end
************************************************************************

      INTEGER OIADTIM,OIADIN,OIREPIT,OIEPSAD,OIEXTIM
      INTEGER SZADTI

      PARAMETER (OIADTIM =  1)
      PARAMETER (OIEXTIM =  2)
      PARAMETER (OIADIN  =  3)
      PARAMETER (OIREPIT =  4)
      PARAMETER (OIEPSAD =  5)
      
      PARAMETER (SZADTI  = 16)

      INTEGER ODTMIN,ODTMAX,ODTFACT,OTIMSTP,OTIMEIN
      INTEGER OEPSADI,OEPSADL,OEPSADU,OEPSNS
      INTEGER SZADTD
      
      PARAMETER (ODTMIN  =  1)
      PARAMETER (ODTMAX  =  2)
      PARAMETER (ODTFACT =  3)
      PARAMETER (OTIMSTP =  4)
      PARAMETER (OTIMEIN =  5)
      PARAMETER (OEPSADI =  6)
      PARAMETER (OEPSADL =  7)
      PARAMETER (OEPSADU =  8)
      PARAMETER (OEPSNS  =  9)

      PARAMETER (SZADTD  = 16)
      