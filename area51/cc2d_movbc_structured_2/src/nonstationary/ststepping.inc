************************************************************************
* This file describes a structure that collects time stepping
* information for the nonstationary Navier Stokes equations.
*
* The structures itselves are realized as double precision
* arrays:
*   DPARAM = array [1..SZTSTD] of double
*          = TTimeStepDParams
* where the elements can be accessed with the help of the offset
* constants defined in this file, e.g. TIMENS=DPARAM(OTIMENS).
*
* The structure defines information about the current simulation
* time, time-stepping information and information about the state
* in the linear extrapolation (if used).
* It is used to configure the behaviour of MGSTP2 how to perform
* one time-step or one substep, respectively.
* For this purpose it contains variables for the different
* weights in each (sub)step of the time-stepping. The structure
* is filled by the main simulation engine NONSTL with data about
* each (sub)step and then passed to MGSTP2 which does the task of
* performing the (sub)step(s).
*
* TTimeStepDParams = record 
*
*   // General information about the current step
*
*   CURTIM : Simulation time of the time step, corresponding to 
*            current solution u_n
*
*   // Theta-scheme
*
*   THETA  : Main configuration parameter of the Theta scheme
*   THETP  : "Adjungated" parameter Theta' of the Theta scheme.
*            Only used in fractional step. 
*            Standard = 0D0 = no fractional step used.
*   ALPHA  : ALPHA-parameter for fractional step. Not used for standard
*            time stepping scheme.
*            Standard = 0D0 = no fractional step used.
*   BETA   : BETA-parameter for fractional step. Not used for standard
*            time stepping scheme.
*            Standard = 0D0 = no fractional step used.
*
*   THSTEP : Theta-scheme identifier for the substeps of the step.
*            =0: Forward Euler, =1: Backward Euler, =0.5: Crank Nicolson
*            Special values for Fractional Step.
*
*   TSTEP  : length of the time step
*
*   // Weights of the differential equation in the substep.
*   // The differential equation in each substep has the form:
*   //
*   //  [ I + TMWGH N(u_n+1) ] u_n+1  +  ... grad (p_n+1) = 
*   //
*   //       [ I + TRMWGH N(u_n+1) ] u_n  +  TR1WGH f_n+1  +  TR2WGH f_n
*   // 
*   // If the simulation has steady (in)homogenuous boundary conditions
*   // and a One-Step Theta-Scheme is used, the RHS build from f_n+1
*   // and f_n can be simplified. In this case:
*   // 
*   //      Theta k f^(n+1) + (1-Theta) k f^n
*   //    = Theta k f^n + (1-Theta) k f^n
*   //    = TRSWGH f^n
*
*   TMWGH  : Weight in th system matrix in this time step
*   TRMWGH : 1st weight in the RHS
*   TR1WGH : 2nd weight in the RHS
*   TR2WGH : 3rd weight in the RHS
*   TRSWGH : alternative RHS weight for steady boundary conditions
*
* end
************************************************************************

      INTEGER OCURTIM
      INTEGER OTHETA,OTHETP,OTHSTEP,OALPHA,OBETA,OTSTEP
      INTEGER OTMWGH,OTRMWGH,OTR1WGH,OTR2WGH,OTRSWGH
      INTEGER OTLEX1,OTLEX2
      INTEGER SZTSTD
      
      PARAMETER (OCURTIM =  1)
      PARAMETER (OTHETA  =  2)
      PARAMETER (OTHETP  =  3)
      PARAMETER (OALPHA  =  4)
      PARAMETER (OBETA   =  5)
      PARAMETER (OTHSTEP =  6)
      PARAMETER (OTSTEP  =  7)
      PARAMETER (OTMWGH  =  8)
      PARAMETER (OTRMWGH =  9)
      PARAMETER (OTR1WGH = 10)
      PARAMETER (OTR2WGH = 11)
      PARAMETER (OTRSWGH = 12)
      PARAMETER (OTLEX1  = 13)
      PARAMETER (OTLEX2  = 14)
      
      PARAMETER (SZTSTD  = 16)
      