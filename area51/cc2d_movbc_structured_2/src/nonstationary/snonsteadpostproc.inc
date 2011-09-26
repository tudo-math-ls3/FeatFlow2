************************************************************************
* This file describes the solver structure for the postprocessing
* routines of the nonstationary Navier-Stokes solver.
*
* It defines two parameter blocks for integer and double parameters,
* respectively. Both parameter blocks contains two parts:
* - An input-block which has to be initialized before the call
*   of the solver
* - A status- & temporary block, which is only valid
*   during the computation and which is maintained by the postprocessing
*   routine internally. They have not to be initialised in advance.
*
* The structures themselves are realized as integer/double precision
* arrays:
*   IPARAM = array [1..SZISPI] of integer 
*          = TInstatPostprocessingIParams
*   DPARAM = array [1..SZISPD] of double
*          = TInstatPostprocessingDParams
*
* The whole structure can be initializes with the INIISP subroutine.
* After initialisation, it's passed to the postprocessing routine of
* the nonstationary solver without being manipulated by the solver;
* it's completely maintained in the postprocessing routine.
*
* TInstatPostprocessingIParams = record 
*
*   // Input parameters: General
*
*   ISAV   : Save intermediate solutions to disc, stepwise.
*             =0: Don't write solution vectors to disc.
*            <>0: Write solution vectors to disc on level |ISAV|.
*                 INSAV>0: Formatted output (standard)
*                 INSAV<0: Unformatted output.
*                 The solutions are written to disc with an 
*                 increasing extension  .1,.2,.3,... until reaching 
*                 .INSVN. Then the counting restarts at 1.
*            What and how solution vectors are written to disc
*            is defined by the output configuration below.
*   INSAV  : Step count for writing solutions to disc in a stepwise
*            manner. This is only active if ISAV<>0. 
*            When INSAV <> 0, the last INSVN accepted solutions 
*            on the finest level are written to disc, regardless of
*            the time stepping.
*            =-1: write everything to disc, predictor steps,
*                 not-accepted corrector steps and accepted corrector
*                 steps.
*            = 0: Write out all accepted solution to disc.
*            > 0: Write out every INSAV'th accepted solution to disc.
*                 Solutions of predictor steps are not saved.
*   INSVN  : Modulo number for how many intermediate solutions are
*            written to disc. When ISAV<>0, INSAV>0, solutions are 
*            written to disc with an increasing extension 
*            .1,.2,.3,... until reaching .INSVN. 
*            Then the counting restarts at 1.
*            The files are written to ./ns/DX.xxxx
*
*   IFINIT : start file number for film output
*
*   IUPSAV : Save intermediate solutions to disc, based on the time-
*            stepping.
*             =0: Don't save time-stepping based intermediate solutions
*                 to disc
*            <>0: Save plain intermediate solutions to disc on level 
*                 |IUPSAV|.
*                 IUPSAV>0: save in formatted form (standard).
*                 IUPSAV<0: save in unformatted form.
*                 The files containing the results are numbered with
*                 increasing extension number .1,.2,.3,...
*                 DUPSAV gives the time interval after which the solution
*                 is to be written out. An equidistant decomposition of
*                 the time scale in time intervals of length DUPSAV
*                 is used, and the solution is written out if the
*                 simulation time reaches such an "output check point".
*                 The files are written to "./film/DUVP.xxxx"
*
*   IFUSAV : Save vertex-based intermediate velocity solutions to disc,
*            based on time-stepping.
*            This works like IUPSAV, but interpolates the solution
*            vectors into vertices and implements boundary conditions
*            into these.
*             =0: Don't write vertex-based intermediate velocity solutions
*                 to disc
*            <>0: Save intermediate velocities to disc on level |IFUSAV|.
*                 IFUSAV>0: Save in formatted form (standard).
*                 IFUSAV<0: Save in unformatted form.
*                 The files are named with increasing extension .1,.2,...
*                 X-velocity is saved to "./film/DU.xxxx"
*                 X-velocity is saved to "./film/DV.xxxx"
*
*   IFPSAV : Save vertex-based intermediate pressure solutions to disc,
*            based on time-stepping.
*            This works like IUPSAV, but interpolates the pressure
*            vectors into vertices and implements boundary conditions
*            into these.
*             =0: Don't write vertex-based intermediate pressure solutions
*                 to disc
*            <>0: Save intermediate pressures to disc on level |IFPSAV|.
*                 IFPSAV>0: Save in formatted form (standard).
*                 IFPSAV<0: Save in unformatted form.
*                 The files are named with increasing extension .1,.2,...
*                 The pressure is saved to "./film/DP.xxxx".
*
*   IFXSAV : Save vertex-based intermediate streamline functions to disc,
*            based on time-stepping.
*            This works like IUPSAV, but calculates the streamfunction
*            and saves this one.
*             =0: Don't write streamfunction to disc
*            <>0: Save streamfunction to disc on level |IFXSAV|.
*                 IFXSAV>0: Save in formatted form (standard).
*                 IFXSAV<0: Save in unformatted form.
*                 The files are named with increasing extension .1,.2,...
*                 The pressure is saved to "./film/DISO.xxxx".
*         
*   IGMV   : Save GMV files to disc, based on time-stepping.
*            =0: Don't save time-stepping based GMV-files to disc.
*            >0: Save GMV-file to disc on level IGMV.
*                The files containing the results are numbered with
*                increasing extension number .1,.2,.3,...
*                DTGMV gives the time interval after which the solution
*                is to be written out. An equidistant decomposition of
*                the time scale in time intervals of length DTGMV
*                is used, and the solution is written out if the
*                simulation time reaches such an "output check point".
*                The files are written to "./gmv/u.gmv.xxxx"
*                with "xxxx" being a counter from IGINIT..infinity
*
*   IGINIT : Start number of file suffix of GMV-files
*
*   ITRCES : Number of tracers per row/rol in each element to generate
*            and simulate.
*
*   IERANA : Cubature rule to use for error analysis in comparison
*            to reference solution.
*            =0: Don't use error analysis
*
*   // Input parameters: Calculation of body forces
*
*   IBDFBD : Cubature formula to use for calculation of body forces
*            (drag/lift) when calculating the forces by line integration
*            on the real boundary.
*            =0: Don't calculate body forces on real boundary
*            >0: Use line cubature formula IBDFBD
*
*   IBDFVI : Cubature formula to use for calculation of body forces
*            (drag/lift) when calculating the forces by volume
*            integration. This is used for the calculation of forces
*            arising on fictitious boundary components.
*            =0: Don't calculate body forces on fictitious boundary obj.
*            >0: Element cubature formula for CB2Q
*            =4: 2x2 Gauss; standard
*
*   IBDFLI : Cubature formula to use for calculation of body forces
*            (drag/lift) when calculating the forces by line
*            integration on the reconstructed interface of fictitious
*            boundary components.
*            =0: Don't calculate body forces on fictitious boundary
*                objects by line integration on the reconstructed
*                interface
*            >0: Element cubature formula for CB1
*            =2: Trapezoidal rule; standard
*
*   // Input parameters. Handles to basic filenames
*   // for file output. If a handle is =0, a standard filename will
*   // be chosen.
*
*   LFLSAV : Handle to basic filename for combined velocity/pressure output,
*            based on number of the time step
*   LFLDUP : Handle to basic filename for combined velocity/pressure output,
*            based on time
*   LFLDU  : Handle to basic filename for film output, X-velocity output
*   LFLDV  : Handle to basic filename for film output, Y-velocity output
*   LFLDP  : Handle to basic filename for film output, pressure output
*   LFLISO : Handle to basic filename for film output, streamline function output
*   LFLGMV : Handle to basic filename of GMV-files
*
*   // Status variables; maintained internally when the solver is 
*   // running; need not to be initialised in advance.
*
*   IFILEN : Number of the solution that is written to disc at last.
*            Used as a counter 1..INSAV for writing solutions
*            to disc every couple of steps.
*
*   IFILN  : Number of the solution that is written to disc at last
*            when writing solutions time-based.
*            Used as a counter 1..infinity for writing solutions
*            to disc every couple of steps.
*   IFILB  : Saves on the beginning of a predictor step a copy
*            of OIFILEN to make it possible to return to a
*            previous time step.
*
*   IGMVN  : Number of the GMV gile that is written to disc at last
*            when writing GMV files time-based.
*            Used as a counter 1..infinity for writing GMV files
*            to disc every couple of steps.
*   IGMVB  : Saves on the beginning of a predictor step a copy
*            of OTGMVN to make it possible to return to a
*            previous time step.
*
*   LTRCS  : Handle to a tracers structure if tracers are
*            maintained for GMV output.
*
*   TRKCNT : Counts how many "tracking" values (drag/lift/...)
*            are written into the tracking files by the standard
*            postprocessing routine.
* end
*
* TInstatPostprocessingDParams = record 
*
*   // Input parameters 
*
*   DUPSAV : The time interval after which solution vectors should be
*            written out; only active if IUPSAV<>0.
*   DTGMV  : The time interval after which GMV files should be
*            written out; only active if IGMV<>0.
*
*   // Status variables; maintained internally when the solver is 
*   // running; need not to be initialised in advance.
*
*   TUPSAV : The point in time when the next solution should be
*            written out to disc if IUPSAV<>0
*   TUPSBK : A backup of TUPSAV, made at the beginning of a corrector
*            step.
*
*   TGMV   : The point in time when the next GMV file should be
*            written out to disc if IGMV<>0
*   TGMVBK : A backup of TGMV, made at the beginning of a corrector
*            step.
* end
************************************************************************

      INTEGER OISAV,OINSAV,OINSVN,OIUPSAV,OIFUSAV,OIFPSAV,OIFXSAV
      INTEGER OIGMV,OIGINIT,OITRCES,OIERANA
      INTEGER OIBDFBD,OIBDFVI,OIBDFLI
      INTEGER OLFLSAV,OLFLDUP,OLFLDU ,OLFLDV ,OLFLDP ,OLFLISO,OLFLGMV
      INTEGER OIFILEN,OIFILN,OIFILB
      INTEGER OIGMVN,OIGMVB,OLTRCS,OTRKCNT,OIFINIT

      INTEGER SZISPI
      
      PARAMETER (OISAV    =  1) 
      PARAMETER (OINSAV   =  2)
      PARAMETER (OINSVN   =  3)
      PARAMETER (OIFINIT  =  4)
      PARAMETER (OIUPSAV  =  5)
      PARAMETER (OIFUSAV  =  6)
      PARAMETER (OIFPSAV  =  7)
      PARAMETER (OIFXSAV  =  8)
      PARAMETER (OIGMV    =  9)
      PARAMETER (OIGINIT  = 10)
      PARAMETER (OITRCES  = 11)
      PARAMETER (OIERANA  = 12) 
      PARAMETER (OIBDFBD  = 13)
      PARAMETER (OIBDFVI  = 14)
      PARAMETER (OIBDFLI  = 15)
      PARAMETER (OLFLSAV  = 16)
      PARAMETER (OLFLDUP  = 17)
      PARAMETER (OLFLDU   = 18)
      PARAMETER (OLFLDV   = 19)
      PARAMETER (OLFLDP   = 20)
      PARAMETER (OLFLISO  = 21)
      PARAMETER (OLFLGMV  = 22)
      PARAMETER (OIFILEN  = 23)
      PARAMETER (OIFILN   = 24)
      PARAMETER (OIFILB   = 25)
      PARAMETER (OIGMVN   = 26)
      PARAMETER (OIGMVB   = 27)
      PARAMETER (OLTRCS   = 28)
      PARAMETER (OTRKCNT  = 29)
      
      PARAMETER (SZISPI   = 48)
      

      INTEGER ODUPSAV,ODTGMV
      INTEGER OTUPSAV,OTUPSBK,OTGMV,OTGMVBK
      
      PARAMETER (ODUPSAV  =  1)
      PARAMETER (ODTGMV   =  2)
      PARAMETER (OTUPSAV  =  3)
      PARAMETER (OTUPSBK  =  4)
      PARAMETER (OTGMV    =  5)
      PARAMETER (OTGMVBK  =  6)
      
      INTEGER SZISPD
      
      PARAMETER (SZISPD   = 16)
      