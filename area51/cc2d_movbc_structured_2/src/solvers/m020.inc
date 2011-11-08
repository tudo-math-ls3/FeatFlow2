************************************************************************
* This file describes the structure of the parameter block of the
* extended multigrid solver M020. It needs basic information about
* the multigrid, therefore the file "CBASICMG.INC" has to be included
* before! Furthermore, as the structure is an extension to the
* general solver structure, the file "SSOLVERS.INC" has to be included
* before.
*
* It defines two parameter blocks for integer and double parameters,
* respectively. Both parameter blocks contain three parts:
* - An input-block which has to be initialized by the caller
* - An output-block which returns the result of the computation
* - A status- & temporary block, which is only valid
*   during the computation and which gives additional information
*   to the callback routines.
* The caller can attach additional elements to the parameter blocks
* with further information for the callback-routines, e.g. matrix-
* handles, flags, tags,...
*
* The structures themselves are realized as integer/double precision
* arrays:
*   IPARAM = array [1..SZ020I] of integer
*          = TM020IParams
*   DPARAM = array [1..SZ020D] of double
*          = TM020DParams
*
* The access to the different variables is organized by
* offset constants defined in this file, e.g.:
*
*   ILEV = IPARAM(OILEV)
*
* These structures are direct extensions to the original solver
* structure, i.e. there has no SZxxxx-constant to be added to
* the offset variable.
* For a proper definition of the constants here, the include files
* CBASICMG.INC and SSOLVERS.INC has to be included prior to this file.
*
* The structures are organized as follows. The first part of the
* structure is as usual a block with general solver parameters.
* Attached to that is a couple of multigrid solver specific
* variables. These only describe the behaviour of the solver
* itself, not of the subcomponents (coarse-grid solver e.g.).
* If such information is necessary, the caller has to attach
* these information to the structure to allow callback-routines
* to access them.
*
* The structure contains some user defined data blocks that allow
* to directly pass information about smoother/grid transfer to
* callback routines. Basically these information blocks provide
* space for 16 integers / doubles each. If the caller needs
* more space, the information has to be attached to the structure
* and the callback routines must be designed properly.
*
* TM020IParams = record
*
*   // The first part of TM020IParams is the whole TSolverIParams-
*   // structure, i.e. TM020IParams is a direct extension of
*   // TSolverIParams:
*
*   STATUS : Result of the solving process.
*            =0: solution computed successfully
*            =1: error in parameters
*            =2: iteration broke down, diverging
*   ITE    : Number of performed multigrid iterations
*   ...
*   NEQ    : Number of equations in the solution/RHS vectors
*
*   // Extended structure for multigrid:
*   // ---------------------------------
*   //
*   // Input parameters, general solution process.
*   // -------------------------------------------
*
*   ICYCLE : Type of multigrid cycle.
*            =0: F-cycle = standard
*            =1: V-cycle
*            =2: W-cycle
*            >2: cycle of higher order
*   NODEFC : Prevent defect checks
*            =0: standard multigrid convergence checks
*            =1: Don't calculate defect norms; always perform
*                NITMAX iterations regardless of the norm of the
*                defect vector.
*                Remark: In this case the results about initial defect,
*                final defect, convergence rate,... are meaningless
*                because they are not calculated.
*   SMTAG : Smoother tag.
*           This (user-defined) tag can be used to inform the
*           callback-routine of the smoother about the type of the
*           actual smoother that should be used.
*           The following smoother tags are typically used:
*           = 0: undefined smoother
*           = 1: Jacobi-smoother
*           = 2: GS-smoother
*           = 3: SOR-smoother
*           = 4: SSOR-smoother
*           = 5: CG-smoother
*           = 6: BiCGStab-smoother
*           = 7: SPAI0-smoother
*           = 8: ILU0-smoother
*           = 9: (M)ILU(s)-smoother
*           =10: SPAI(s)-smoother
*           =11: Multigrid-smoother
*           >50: user-defined smoother tag
*           This tag has can be initialized by the initialization
*           routines and be used in the callback-routines to
*           choose a smoother. The multigrid algorithm itself does
*           not use the information stored here.
*
*           Identifiers 1..4 are designed for standard smoother
*           routines, these routines are configured only by one
*           parameter. In the TM020DParams structure there is
*           a double precision tag OMGSM reserved to configure
*           such a "primitive" smoother.
*           As long as less than 16 integers and 16 doubles are
*           necessary to configure a smoother, the initialization/
*           callback routine can also use the SMTHI/SMTHD blocks
*           in the multigrid structure to store configuration
*           parameters for the smoother.
*           All other smoothers configured by more than one
*           piece of information (e.g. matrix handles,...) need
*           more space to store their configuration. In this
*           case additional information has to be added to this
*           multigrid structure.
*           In all cases the callback routine must be written
*           appropriately to fetch the information from the MG
*           structure to correctly perform its task.
*
*   SMTHI : array [1..16] of integer
*            Additional parameters to configure the smoother.
*            User-defined, not used by the multigrid algorithm
*            internally.
*            As a convention, if the number of smoothing steps is
*            constant on each level:
*            - if symmetrical smoothing is used, SMTHI(1) should be
*              initialized to the number of smoothing steps
*            - if the number of pre- and postsmoothing steps differ,
*              SMTHI(1) should be initialized to the number of pre-
*              and SMTHI(2) to the number of postsmoothing steps.
*
*   PRRSI : array [1..16] of integer
*            Additional parameters to configure the grid transfer.
*            User-defined, not used by the multigrid algorithm
*            internally.
*
*   UPARI : array [1..16] of integer
*            Additional general parameters block.
*            User-defined, not used by the multigrid algorithm
*            internally.
*
*   // Output-Tags for messages on the screen. MSGTRM is the output
*   // level for the multigrid. The other output levels are ignored
*   // by the MG-algorithm but should be used by the callback routines
*   // to determine their level of output.
*
*   MSGTRM : Level of output to terminal of MG-algorithm.
*            =0: no output
*            =1: basic output
*            =2: standard output = standard
*   MSMTRM : Level of output to terminal for smoother
*            =0: no output = standard
*            =1: basic output
*            =2: extended output
*   MSTTRM : Level of output to terminal for step length control
*            =0: no output = standard
*            =1: basic output
*            =2: extended output
*   MCGTRM : Level of output to terminal for coarse grid solver
*            =0: no output = standard
*            =1: basic output
*            =2: extended output
*   MBCTRM : Level of output to terminal for boundary implementation
*            =0: no output = standard
*            =1: basic output
*            =2: extended output
*
*   // Input-parameters, problem-specific:
*   // (depending on the actual situation)
*
*   NLMIN  : minimum level, where the coarse grid solver should solve;
*            >= 1.
*   NLMAX  : maximum level, where the solution should be computed;
*            <= NNLEV!
*   KNEQ   : array [1..NNLEV] of integer
*            Number of equations on each level.
*   KOFFX  : array [1..NNLEV] of integer
*            Starting indices of the solution vectors in DX() on the
*            different levels
*   KOFFB  : array [1..NNLEV] of integer
*            Starting indices of the RHS vectors in DX() on the
*            different levels
*   KOFFD  : array [1..NNLEV] of integer
*            Starting indices of temporary vectors in DX() on the
*            different levels
*   KPRSM  : array [1..NNLEV] of integer
*            Number of presmoothing steps on each level
*   KPOSM  : array [1..NNLEV] of integer
*            Number of postsmoothing steps on each level
*
*   // Status variables, only valid during the solution process
*   // --------------------------------------------------------
*
*   ILEV   : Current level in multigrid
*
*   KCBAUX : Starting address of a vector array [1..NEQ] of double
*            in DD, NEQ=NEQ(finest grid). (More precisely: DD(NLMAX).)
*            This array represents an auxiliary memory
*            block DX/DB/DD(KCBAUX:KCBAUX+NEQ-1) that
*            can be used by the callback routines DPROL, DSTEP
*            and DEX for intermediate calculations.
*            The restriction routine DREST can only use it if its fine
*            grid is < NLMAX or after the restriction is completed.
*            DEXS is not allowed to use it. The other routines (DFILT,
*            DAX,...) can use it only if the current level is < NLMAX.
*            The content of this array is undefined on call of the
*            callback routine and will be treated as undefined by MG
*            when the callback routine leaves.
*            =0 : KBCAUX not available.
*            <>0: KBCAUX available.
*
*   KIT    : array [1..NNLEV] of integer
*            Internal temporary variables.
*   KIT0   : array [1..NNLEV] of integer
*            Internal temporary variables.
* end
*
* TM020DParams = record
*
*   // The first part of TM020DParams is the whole TSolverIParams-
*   // structure, i.e. TM020DParams is a direct extension of
*   // TSolverPParams:
*
*   DEFINI : initial defect
*   DEFFIN : final defect after last iteration
*   ...
*   PRECD  : array [1..16] of integer
*
*   // Output variables
*   // ---------------------------------
*   // There are no additional output variables added to the standard
*   // output variables of the general solver structure
*
*   // Extended structure for multigrid:
*   // ---------------------------------
*   // Timing output parameters; only set if IMTIM > 0:
*
*   TMMG   : time necessary for the multigrid algorithm internally
*   TMPROL : time for prolongation
*   TMREST : time for restriction
*   TMDEF  : time for defect calculation
*   TMSMTH : time for smoothing operations
*   TMCGSL : time for coarse-grid solver
*   TMBC   : time for implementation of boundary conditions
*   TMCGC  : time for coarse-grid corrections
*
*   // Input parameters
*   // ---------------------------------
*
*   STPMIN : Mimimum step length for optimal coarse grid correction.
*            Not used by the MG solver internally, but can be used
*            by the callback routines for step length control.
*            Standard = 1D0
*   STPMAX : Maximum step length for optimal coarse grid correction.
*            Not used by the MG solver internally, but can be used
*            by the callback routines for step length control.
*            Standard = 1D0
*   OMGSM  : user-defined parameter to configure the smoother
*   SMTHD  : array [1..16] of double
*            Additional parameter to configure the smoother.
*            User-defined.
*   PRRSD  : array [1..16] of double
*            Additional parameters to configure the grid transfer.
*            User-defined, not used by the multigrid algorithm
*            internally.
*   UPARD  : array [1..16] of double
*            Additional general parameters block.
*            User-defined, not used by the multigrid algorithm
*            internally.
* end
*
* The standard settings of the parameters can be initialised with
* INM020. This routine calls INGSLV to initialize the general solver
* parameters, too.
************************************************************************

***************
* TM020IParams
***************

      INTEGER ONLMIN ,ONLMAX,OICYCLE
      INTEGER OKNEQ ,OKOFFX,OKOFFB,OKOFFD
      INTEGER OILEV ,OKIT  ,OKIT0, OKPRSM, OKPOSM, ONODEFC
      INTEGER OMSMTRM,OMSTTRM,OMCGTRM,OMBCTRM,OSMTAG,OSMTHI,OPRRSI
      INTEGER OUPARI,OKCBAUX
      INTEGER SZ020I

C The constants 1..SZSLVI are defined in SSOLVERS.INC.
C The structure is a direct extension to the TSolverIParams structure,
C so the first variable starts right after the SZSLVI'th element.
      
      PARAMETER (OICYCLE = SZSLVI+ 1)
      PARAMETER (ONODEFC = SZSLVI+ 2)
      PARAMETER (OSMTAG  = SZSLVI+ 3)
      PARAMETER (OSMTHI  = SZSLVI+ 4)
      PARAMETER (OPRRSI  = SZSLVI+ 20)
      PARAMETER (OUPARI  = SZSLVI+ 36)
      PARAMETER (OMSMTRM = SZSLVI+ 52)
      PARAMETER (OMSTTRM = SZSLVI+ 53)
      PARAMETER (OMCGTRM = SZSLVI+ 54)
      PARAMETER (OMBCTRM = SZSLVI+ 55)
      
      PARAMETER (ONLMIN  = SZSLVI+ 56)
      PARAMETER (ONLMAX  = SZSLVI+ 57)
      PARAMETER (OKNEQ   = SZSLVI+ 58)
      PARAMETER (OKOFFX  = OKNEQ+NNLEV)
      PARAMETER (OKOFFB  = OKOFFX+NNLEV)
      PARAMETER (OKOFFD  = OKOFFB+NNLEV)
      PARAMETER (OKPRSM  = OKOFFD+NNLEV)
      PARAMETER (OKPOSM  = OKPRSM+NNLEV)
      
      PARAMETER (OILEV   = OKPOSM+NNLEV)
      PARAMETER (OKCBAUX = OILEV+1)
      PARAMETER (OKIT    = OKCBAUX+1)
      PARAMETER (OKIT0   = OKIT+NNLEV)
      
      PARAMETER (SZ020I  = OKIT0+NNLEV+15)
                 
C OKIT0..SZ020I : reserved
C SZ020I+1..*   : user defineable
                 
***************
* TM020DParams
***************
                 
      INTEGER OTMMG  ,OTMPROL,OSTPMIN,OSTPMAX
      INTEGER OTMREST,OTMDEF ,OTMSMTH,OTMCGSL,OTMBC,OTMCGC,OSMTHD,OOMGSM
      INTEGER OPRRSD ,OUPARD
      INTEGER SZ020D
                 
C The constants 1..SZSLVD are defined in SSOLVERS.INC.
C The structure is a direct extension to the TSolverDParams structure,
C so the first variable starts right after the SZSLVD'th element.

      PARAMETER (OTMMG   = SZSLVD+ 1)
      PARAMETER (OTMPROL = SZSLVD+ 2)
      PARAMETER (OTMREST = SZSLVD+ 3)
      PARAMETER (OTMDEF  = SZSLVD+ 4)
      PARAMETER (OTMSMTH = SZSLVD+ 5)
      PARAMETER (OTMCGSL = SZSLVD+ 6)
      PARAMETER (OTMBC   = SZSLVD+ 7)
      PARAMETER (OTMCGC  = SZSLVD+ 8)

      PARAMETER (OSTPMIN = SZSLVD+ 9)
      PARAMETER (OSTPMAX = SZSLVD+10)

      PARAMETER (OOMGSM  = SZSLVD+11)
      PARAMETER (OSMTHD  = SZSLVD+12)
      
      PARAMETER (OPRRSD  = SZSLVD+28)
      PARAMETER (OUPARD  = SZSLVD+44)
      
C ODSMPAR+16..SZ020D: reserved
      
      PARAMETER (SZ020D  = OSMTHD+64)

C SZ020D+1..* : user-defineable
