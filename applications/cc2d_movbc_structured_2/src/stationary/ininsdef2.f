************************************************************************
* This file collects the initialization routines for NSDEF.
*
* The routines here build the general parameter structures for the
* call to NSDEF2 solver.
************************************************************************

************************************************************************
* Initialize NSDEF stationary solver structures
*
* This initializes the IPARAM/DPARAM structure arrays, which are
* used in the call of NSDEF, with standard values. If desired,
* values of the COMMON-blocks are transferred to these structure
* arrays.
*
* The routine does not allocate any memory, it just initializes the
* parameters.
*
* In:
*   IPARAM : array [1..SZNSDI] of integer
*            Integer parameter block for NSDEF
*   DPARAM : array [1..SZNSDD] of double
*            Double precision parameter block for NSDEF
*   IC2PAR : =0: initialize IPARAM/DPARAM only with standard values
*            =1: initialize IPARAM/DPARAM with standard values and
*                transfer values of COMMON-blocks into them
*
* Out:
*   IPARAM,
*   DPARAM : Initialized structures
*
* This routine does not allocate any memory and does not initialize
* the DString-handles of the basic filenames for file output. This
* has to be done by the caller if necessary!
************************************************************************

      SUBROUTINE INISTS(IPARAM,DPARAM,IC2PAR)
      
      IMPLICIT NONE
      
      INCLUDE 'cout.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'stria.inc'
      INCLUDE 'ssolvers.inc'
      INCLUDE 'slinsol.inc'
      INCLUDE 'slinsol2dns.inc'
      INCLUDE 'smat2dns.inc'
      INCLUDE 'stiming.inc'
      INCLUDE 'snsdef.inc'
      
      INCLUDE 'cnonlinsol_cc2d.inc'
      INCLUDE 'cpostproc.inc'
      
      INTEGER IPARAM(SZNSDI),IC2PAR
      DOUBLE PRECISION DPARAM(SZNSDD)
      
C     Clean the arrays
      
      CALL LCL1 (DPARAM,SZNSDD)
      CALL LCL3 (IPARAM,SZNSDI)
      
C     Set up standard values.
C     At first call the standard initialization routine for the standard
C     solver structure.

      CALL INGSLV (IPARAM, DPARAM)

C     Some standard values are not used by NSDEF, set them to 0.
C     Others have different <> 0 standard values - initialize them.

      IPARAM (ONITMAX) = 20
      IPARAM (OTIM  )  = 2
      IPARAM (OIASRLN) = 0
      IPARAM (OMSGTRM) = 2
      IPARAM (OMTRMRS) = 1

      DPARAM (OEPSREL) = 0D0
      DPARAM (OEPSABS) = 0D0
      DPARAM (ODIVREL) = 0D0
      DPARAM (ODIVABS) = 0D0
      DPARAM (OVECZER) = 0D0
      
C     Initialize extended parameters which are different from 0.

      IPARAM (OSLTAG)  = 20
      IPARAM (OMSGLIN) = 1
      
C     Should we transfer COMMON-block variables?

      IF (IC2PAR.NE.0) THEN
      
        IPARAM (ONITMIN) = INLMIN
        IPARAM (ONITMAX) = INLMAX
        IPARAM (OMSGTRM) = MT
        
        IPARAM (OIAUSAV) = IAUSAV
        
        DPARAM (OOMGINI) = OMGINI
        DPARAM (OOMGMIN) = OMGMIN
        DPARAM (OOMGMAX) = OMGMAX
              
C       Transfer stopping criteria
      
        DPARAM (OEPSD  ) = EPSD  
        DPARAM (OEPSDIV) = EPSDIV
        DPARAM (OEPSUR ) = EPSUR 
        DPARAM (OEPSPR ) = EPSPR 
        DPARAM (ODMPD  ) = DMPD  
      
      END IF
      
      END
      
************************************************************************
* Initialize structure for linear solver in NSDEF stationary solver
*
* This initializes the IPARAM/DPARAM structure arrays, which define the
* behaviour of the linear solver that is used as sub-solver in
* NSDEF2. If desired, values of the COMMON-blocks are transferred to 
* these structure arrays.
*
* The routine does not allocate any memory, it just initializes the
* parameters.
*
* IPARAM/DPARAM should be filled with 0 before calling this routine,
* as here only the nonzero parameters are set.
*
* In:
*   NLMIN  : Minimum level in multigrid; >= 1
*   NLMAX  : Maximum level in multigrid; <= 9
*   IPARAM : array [1.SZ020I+3*SZSLVI] of integer
*            Integer parameter block of solver structure for NSDEF2.
*   DPARAM : array [1..SZ020D+3*SZSLVD] of double
*            Double precision parameter block of structure 
*            of linear solver
*   ISLV   : Type of solver that should be used.
*            =0: Multigrid solver
*            =1: UMFPACK (1-level) solver
*   IC2PAR : =0: initialize IPARAM/DPARAM only with standard values
*            =1: initialize IPARAM/DPARAM with standard values and
*                transfer values of COMMON-blocks into them
*
* Out:
*   IPARAM,
*   DPARAM : Initialized structures for solver, multigrid, smoother and
*            coarse grid solver
*
* After calling this routine, the first solver structure in IPARAM
* decides on which solver is used. The caller should call the
* correct solver depending on the solver tag in this structure.
************************************************************************

      SUBROUTINE INLSST (NLMIN,NLMAX,IPARAM,DPARAM,ISLV,IC2PAR)
      
      IMPLICIT NONE
      
      INCLUDE 'cout.inc'
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'ssolvers.inc'
      INCLUDE 'm020.inc'
      
      INCLUDE 'clinsol_cc2d.inc'
      
C     parameters
      
      INTEGER IPARAM(*),IC2PAR,NLMIN,NLMAX,ISLV
      DOUBLE PRECISION DPARAM(*)
      
C     Should we use a standard MG solver?

      IF (ISLV.EQ.0) THEN
      
C       Ok. Call INMGST to initialize the standard MG solver
C       structure. IPARAM contains two or three substructures:
C       - Multigrid solver with VANCA smoother
C       - (Vanca Coarse) grid solver
C       or (in future versions)
C       - Multigrid solver with VANCA smoother
C       - Smoother-structure for advanced
C       - (Vanca Coarse) grid solver

        CALL INMGST (NLMIN,NLMAX,IPARAM,DPARAM,IC2PAR)
      
      ELSE IF (ISLV.EQ.1) THEN
      
C       The user wants to use BiCGStab-solver with a
C       multigrid preconditioner.
C       IPARAM contains three or four substructures:
C       - BiCGStab-solver
C       - Multigrid preconditioner with VANCA smoother
C       - (Vanca Coarse) grid solver
C       or (in future versions)
C       - BiCGStab-solver
C       - Multigrid preconditioner 
C       - Smoother-structure for advanced
C       - (Vanca Coarse) grid solver
C
C       Initialize the BiCGStab-structure with standard
C       solver parameters:

        CALL INGSLV (IPARAM, DPARAM)
        
C       Initialize the tag to indicate a BiCGStab-solver:

        IPARAM(OSLTAG) = 6

C       Call INMGST to initialize the standard MG solver
C       structure, which is later modified to be a 
C       preconditioning structure. The structure is found
C       at the second block of IPARAM/DPARAM behind the
C       BiCGStab solver structure.

        CALL INMGST (NLMIN,NLMAX,
     *               IPARAM(SZSLVI+1),DPARAM(SZSLVD+1),IC2PAR)
     
C       Modify the solver structure to form the structure of a
C       multigrid preconditioner.
C       Initialize the tag to indicate a Multigrid preconditioner.
C       Inform the BiCGStab solver that we have a preconditioner.

        IPARAM(SZSLVI+OSLTAG) = 111
        IPARAM(OPCTAG)        =  11
        
C       Initialize minimum and maximum iterations of the 
C       preconditioner:
      
        IPARAM(SZSLVI+ONITMIN) = 2
        IPARAM(SZSLVI+ONITMAX) = 2

C       Reduce the message level of the MG components by 2

        IPARAM(SZSLVI+OMSGTRM) = IPARAM(SZSLVI+OMSGTRM)-2
        IPARAM(SZSLVI+OMCGTRM) = IPARAM(SZSLVI+OMCGTRM)-2
        IPARAM(SZSLVI+OMSMTRM) = IPARAM(SZSLVI+OMSMTRM)-2
        
C       Set the message level for BiCGStab differently to the
C       standard initialization.

        IPARAM(OMSGTRM) = MT-2

C       Transfer COMMON block variables from DAT file
        
        IF (IC2PAR.NE.0) THEN
        
C         Number of preconditioning steps is written to min.+max. 
C         number of iterations of the MG component.
C         This deactivates any stopping criterion.
        
          IPARAM (SZSLVI+ONITMIN) = IPCMG
          IPARAM (SZSLVI+ONITMAX) = IPCMG
          
C         Min./Max. number of iterations of BiCGStab
          
          IPARAM (ONITMIN) = ILMIN
          IPARAM (ONITMAX) = ILMAX

C         Stopping criteria for BiCGStab.

          DPARAM (OEPSREL) = EPSMG
          DPARAM (OEPSABS) = DMPMG
        END IF
        
      ELSE IF (ISLV.EQ.2) THEN
      
C       The user wants to use a standard UMFPACK4 direct solver.
C       We can easily prepare the structure in IPARAM, which
C       contains only one block, configured as standard solver
C       structure:

        CALL INGSLV (IPARAM, DPARAM)
        
C       Note the type of solver in IPARAM

        IPARAM(OSLTAG) = 13
        
C       Nothing more to do.
      
      END IF

      END

************************************************************************
* Initialize Multigrid structure for NSDEF-solver
*
* This initializes the IPARAM/DPARAM structure arrays, which define the
* behaviour of the multigrid solver that is used as sub-solver in
* NSDEF2. If desired, values of the COMMON-blocks are transferred to 
* these structure arrays.
*
* The routine does not allocate any memory, it just initializes the
* parameters.
*
* IPARAM/DPARAM should be filled with 0 before calling this routine,
* as here only the nonzero parameters are set.
*
* In:
*   NLMIN  : Minimum level in multigrid; >= 1
*   NLMAX  : Maximum level in multigrid; <= 9
*   IPARAM : array [1.SZ020I+3*SZSLVI] of integer
*            Integer parameter block of MG structure for NSDEF2
*   DPARAM : array [1..SZ020D+3*SZSLVD] of double
*            Double precision parameter block of MG structure 
*            for NSDEF2
*   IC2PAR : =0: initialize IPARAM/DPARAM only with standard values
*            =1: initialize IPARAM/DPARAM with standard values and
*                transfer values of COMMON-blocks into them
*
* Out:
*   IPARAM,
*   DPARAM : Initialized structures for multigrid, smoother and
*            coarse grid solver
************************************************************************

      SUBROUTINE INMGST (NLMIN,NLMAX,IPARAM,DPARAM,IC2PAR)
      
      IMPLICIT NONE
      
      INCLUDE 'cout.inc'
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'ssolvers.inc'
      INCLUDE 'm020.inc'
      
      INCLUDE 'clinsol_cc2d.inc'
      
C     parameters
      
      INTEGER IPARAM(*),IC2PAR,NLMIN,NLMAX
      DOUBLE PRECISION DPARAM(*)

C     local variables
      
      INTEGER IDX
      
C     The xPARAM structure arrays contain information for three
C     different type of iterative algorithms: (Multigrid) solver,
C     smoother and coarse grid solver.
C
C     The standard configuration is:
C       - Multigrid solver
C       - Vanca Smoother
C       - Vanca Coarse grid solver
C     So we initialize the structures according to that.
C
C     The first part contains information that configures the 
C     solver itself. Initialize it with the standard MG initialization
C     routine:

      CALL INM020 (IPARAM, DPARAM)
      
C     The second part configures the smoother. Initialize that with
C     standard values:

      CALL INGSLV (IPARAM(1+SZ020I), DPARAM(1+SZ020D))
      
C     The third part configures the coarse grid solver, which is in
C     the standard case a VANCA pseudo-solver. Nevertheless, it's an 
C     iterative  solver that requires an appropriate solver structure.

      CALL INGSLV (IPARAM(1+SZ020I+SZSLVI), DPARAM(1+SZ020D+SZSLVD))
      
      IF (IC2PAR.NE.0) THEN
      
C       The user wants us to transfer COMMON block variables to
C       the solver structure for initialization. This allowes to
C       configure/initialize the solver more precisely.

C       Message level is two below the current message level

        IPARAM(OMSGTRM) = MT-2

C       Set the output level of the coarse grid solver to the 
C       general output level - 4

        IPARAM(OMCGTRM) = MT-4

C       Set the output level of the smoother to the 
C       general output level - 6

        IPARAM(OMSMTRM) = MT-6

C       Preset of minimum/maximum level. Might be changed in NSDEF.

        IPARAM(ONITMIN) = ILMIN
        IPARAM(ONITMAX) = ILMAX
        
C       Stopping criteria for MG solver:

        DPARAM(OEPSABS) = EPSMG
        DPARAM(OEPSREL) = DMPMG
        
C       Minimum/maximum value for the optimal correction:

        DPARAM(OSTPMIN) = AMINMG
        DPARAM(OSTPMAX) = AMAXMG

C       Configure the smoother. We only support basic smoothing,
C       so we don't have to pass much information to the callback 
C       routines. At first, initialize the smoother tag to identify
C       out VANCA. That's not used by the MG algorithm internally,
C       but it's good style to do this...
C
C       VANCA = 50+ISM

        IPARAM(OSMTAG) = 50+ISM
        
C       Initialize the relaxation parameter of the smoother

        DPARAM(OOMGSM) = RLXSM        

C       Build up the number of smoothing steps using the parameters
C       from the COMMON blocks.

        DO IDX = NLMIN,NLMAX
          IPARAM (OKPRSM+IDX-1) = NSM * NSMFAC**(NLMAX-IDX)
          IPARAM (OKPOSM+IDX-1) = NSM * NSMFAC**(NLMAX-IDX)
        END DO
                  
C       That's it for the smoother. Now configure the grid transfer.
C       We configure PRRSI/PRRSD as follows to pass parameters to the
C       callback routines for prolongation/restriction:
C
C         PRRSI(1) = IINT
C         PRRSI(2) = IAPRM
C         PRRSI(3) = IAVPR
C
C         PRRSD(1) = DPREP

        IPARAM(OPRRSI)   = IINT
        IPARAM(OPRRSI+1) = IAPR
        IPARAM(OPRRSI+2) = IAVPR

        DPARAM(OPRRSD)   = DPREP
        
C       Now we come to the initialization of the coarse grid solver.
C
C       Remember that all variables initialized here for the smoother/
C       coarse grid solver are only pseudotypically prescribed!
C       They are initialized in a rather structured way, but only serve
C       to pass information to the corresponding callback-routines of NSDEF!
C       All information written in the structures here is only used
C       by the callback-routines of NSDEF - nowhere else!
C       That way we have a template so that we can theoretically switch
C       to another smoother/solver later, as the structure is the same
C       for all kinds of smoothers/solvers.

C       Configure the coarse grid solver:

        IDX = SZ020I
        
C       Maximum number of iterations:

        IPARAM(IDX+ONITMAX) = NSL
        
C       Which type of coarse grid solver do we actually have? 
C       ISL decides about that. 

        IF (ISL.EQ.10) THEN
        
C         ISL=10 -> Type of solver: 13=UMFPACK4

          IPARAM(IDX+OSLTAG) = 13
          
        ELSE 
                
C         ISL=0 -> Type of solver: 50=VANCA
C         ISL=1 -> Type of solver: 51=Extended VANCA
C         ISL=2 -> Type of solver: 50=VANCA, defect correction
C         ISL=3 -> Type of solver: 51=Extended VANCA, defect correction
          
          IPARAM(IDX+OSLTAG) = 50+ISL
        
        END IF
      
C       Absolute and relative error bound

        IDX = SZ020D
        
        DPARAM(IDX+OEPSABS) = EPSSL
        DPARAM(IDX+OEPSREL) = DMPSL
        
C       Standard relative divergence criterion

        DPARAM(IDX+ODIVREL) = 1D6
                
C       Save the OMEGA relaxation parameter for coarse grid solver
        
        DPARAM(IDX+OOMEGA)  = RLXSL
 
      END IF
      
C     Activate the following setting to get coarse grid solver output:
C       
C      IPARAM(OMCGTRM) = 2
        
      END
