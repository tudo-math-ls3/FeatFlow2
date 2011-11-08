************************************************************************
* This file describes the structures of the standard parameter block
* for solvers, extended smoothers and extended preconditioners.
* This structures has to be extended to a specialized
* structure, which contains the solver specific parameters. There are
* two structures, an integer and a double precision structure.
*
* Both parameter blocks contain three parts:
* - An input-block which has to be initialized by the caller
* - An output-block which returns the result of the computation
* - A status- & temporary block, which is only valid
*   during the computation and which gives additional information
*   to the callback routines.
*
* The caller can attach additional elements to the parameter blocks
* with further information for the callback-routines, e.g. matrix-
* handles, flags, tags,... There are intentionally no "exchange"
* parameters inside of the structure to prevent unpredictable
* results (like in former solver calls, where the ITE-parameter
* was on call the minimum number of iterations and returned the
* actual performed iteration count) when calling the routine
* more than once. This allows to call the solver more
* than one time with the same structure, getting the same results.
*
* The structures themselves are realized as integer/double precision
* arrays:
*   IPARAM = array [1..SZSLVI] of integer
*          = TSolverIParams
*   DPARAM = array [1..SZSLVD] of double
*          = TSolverIDParams
*
* The access to the different variables is organized by
* offset constants defined in this file, e.g.:
*
*   STATUS = IPARAM(OSTATUS)
*
* The usage of this solver structure in combination with different
* solvers is typically the following:
* - Take the general solver structure and extend it to a solver-
*   specific structure by adding further solver-specific variables
*   to it for input/output/status
* - Initialise the specialised solver structure by a solver-specific
*   initialisation routine (INxxxx). Initialise the general solver
*   structure (e.g. by INGSLV) - this is normally done by the
*   specialized initialization routine, too.
* - Call the solver with the structures as parameters
* - When the solver is finished, extract the results of the solver
*   from the output-blocks of the specialized solver structure
* This means, for a standard solver, the specialized solver structure
* is designed to have the following form:
*   record
*     - General output variables
*     - General input parameters
*     - General status variables
*     - Specific output variables
*     - Specific input parameters
*     - Specific status variables
*   end
* The access to these variables is done with the constants as described
* above. The SZxxxx-constant describes the size of the structure.
* Standard solver parameters contain in each structure a set
* of 99 elements. The SZxxxx-constants are designed to cover 128
* elements. This allowes solvers with a limited set of additional
* parameters to include these directly in the array without
* defining a new SZxxxx-constant for that solver. The additional
* data can be stored directly at position 100..128 in the array.
* Solvers with more parameters need to attach complete blocks
* and have tro define a new SZxxxx-constant for that block.
*
* To provide direct support for standard matrix-vector multiplication
* routines or standard preconditioning routines are used, the caller
* can provide the appropriate information in the MATRx and PRECx
* subarrays of the structure. The standard matrix-vector multiplication
* and preconditioning routines are designed to access these blocks
* inside if the IPARAM/DPARAM block directly. For more complex
* matrix/vector multiplication or preconditioning, the caller has to
* attach the information to the end of the parameter block and
* create appropriate callback-routines to access this information.
*
* If some solvers are nested (like e.g. using the multigrid
* solver), the specific variable blocks in the total structure might
* even contain variables for the sub-solvers. As an example we describe
* here an example of a possible structure of the multigrid-specific
* solver structure:
*   record
*     - General output variables
*     - General input parameters
*     - General status variables
*     - Specific output variables, only valid for MG-solver
*     - Specific input parameters, only valid for MG-solver
*     - Specific status variables, only valid for MG-solver
*
*     - General output variables for coarse-grid solver
*     - General input parameters for coarse-grid solver
*     - General status variables for coarse-grid solver
*     - Specific output variables for coarse-grid solver
*     - Specific input parameters for coarse-grid solver
*     - Specific status variables for coarse-grid solver
*   end
* As it can be seen, the total block for the multigrid solver decomposes
* into two blocks. The first general/specific block describes the
* parameters for the multigrid solver itself. Attached to that there
* is a second general/specific block that contains the parameters
* for the coarse grid solver:
* As the coarse grid solver is itself derived from a standard solver,
* it needs to be configured by such a block. This block is now embedded
* into the multigrid structure. This allows the coarse-grid solver
* callback-routine inside of the multigrid solver to directly
* call the solver component with the appropriate block without building
* it every time, which is the typical usage.
* It's even possible to attach more of such solver-blocks to the
* multigrid block in case that this is necessary (e.g. more than one
* coarse grid solver, and the callback-routine switches dynamically
* between them), or to attach blocks of information about computational
* grids to it. To access the different structure blocks,
* the SZxxxx-size constants have to be used to calculate the
* correct offset of the block inside of the whole structure!
* E.g. in the above example if the M020-solver is used, the parameters
* for the coarse grid solver are located at IPARAM(SZM020I+1) and
* DPARAM(SZM020I+1), respectively.
*
* The basic structures are organized as below. The specific
* solver might ignore some of the variables if necessary. In this
* case these variables have to be set to 0 by the initialization
* routine and haven't to be tackled by the solver (except where
* mentioned).
*
* TSolverIParams = record
*
*   // Output parameters
*
*   STATUS : Result of the solving process.
*            =0: solution computed successfully
*            =1: error in parameters
*            =2: iteration broke down, diverging
*
*   ITE    : Number of performed iterations, if the solver
*            is of iterative nature.
*            Must be set to 1 by the solver if not used
*            (indicating at least 1 performed iteration, which is
*             always the case).
*
*   // Input parameters
*
*   SLTAG  : Solver tag
*            This us a user-defineable tag which can be used
*            to identify the solver/smoother/preconditioner, the structure
*            is set up to. Should be set by the initialization routine
*            of this structure to a unique identifier to allow callback
*            routines to decide, for which specific component this structure
*            is prepared for.
*            Can be set to 0 if not used.
*            The following identifiers are reserved:
*            =   0: unknown solver structure
*            =   1: Jacobi-solver
*            =   2: GS-solver
*            =   3: SOR-solver
*            =   4: SSOR-solver
*            =   5: CG-solver
*            =   6: BiCGStab-solver
*            =   7: SPAI0-solver
*            =   8: ILU0-solver
*            =   9: (M)ILU(s)-solver
*            =  10: SPAI(s)-solver
*            =  11: Multigrid-solver
*            =  12: UMFPACK2-solver
*            =  13: UMFPACK4-solver
*
*            = 100: unknown preconditioner
*            = 101: Jacobi-preconditioner
*            = 102: GS-preconditioner
*            = 103: SOR-preconditioner
*            = 104: SSOR-preconditioner
*            = 105: CG-preconditioner
*            = 106: BiCGStab-preconditioner
*            = 107: SPAI0-preconditioner
*            = 108: ILU0-preconditioner
*            = 109: (M)ILU(s)-preconditioner
*            = 110: SPAI(s)-preconditioner
*            = 111: Multigrid-preconditioner
*
*            = 200: unknown smoother
*            = 201: Jacobi-smoother
*            = 202: GS-smoother
*            = 203: SOR-smoother
*            = 204: SSOR-smoother
*            = 205: CG-smoother
*            = 206: BiCGStab-smoother
*            = 207: SPAI0-smoother
*            = 208: ILU0-smoother
*            = 209: (M)ILU(s)-smoother
*            = 210: SPAI(s)-smoother
*            = 211: Multigrid-smoother
*
*            As the general solver structure does not only serve as a basis
*            for solvers but also as a basis for smoothers and preconditioners,
*            the reserved tag numbers separate into different ranges. All
*            numbers > 1000 can be defined by the user of this structure
*            to define other components as the predefined ones.
*            (Remember: This tag is only an info for the callback routines
*             to decide how to react! It must be set and used problem-
*             specifically, so the above numbers are only suggestions!)
*   NITMIN : minimum number of iterations, if the solver is of
*            iterative nature; otherwise not used.
*            >= 0.
*   NITMAX : maximum number of iterations, if the solver is of
*            iterative nature; otherwise not used.
*            >= 0.
*   IFILT  : Filtering tag.
*            This (user-defined) tag can be used to inform the
*            callback-routine for the filtering about the type of
*            the actual filter that should be used.
*            Also informs the solver about whether the callback-routine
*            of the filter should be called at all (if the solver
*            supports filtering).
*            Must be set to 0 if not used. If <> 0, the callback-routine
*            for the filter is called to perform the filtering.
*            =0: don't use filtering; subroutine DFILT is not called
*            =1: perform filtering: frequently call DFILT
*            The filtering routine DFILT itself is solver-specific;
*            see the documentation of the appropriate solver about
*            the syntax of the filtering routine for that solver.
*   OPCTAG : Preconditioner tag.
*            This (user-defined) tag can be used to inform the
*            callback-routine of the preconditioner about the type of
*            the actual preconditioner that should be used.
*            Also informs the solver about whether the callback-routine
*            of the preconditioner should be called at all
*            (if the solver supports preconditioning).
*            Must be set to 0 if not used. If <> 0, the callback-routine
*            for the preconditioner is called to perform the
*            preconditioning.
*            = 0: no preconditioning, don't call preconditioner
*                 callback-routine
*            = 1: Jacobi-preconditioner
*            = 2: GS-preconditioner
*            = 3: SOR-preconditioner
*            = 4: SSOR-preconditioner
*            = 5: CG-preconditioner
*            = 6: BiCGStab-preconditioner
*            = 7: SPAI0-preconditioner
*            = 8: ILU0-preconditioner
*            = 9: (M)ILU(s)-preconditioner
*            =10: SPAI(s)-preconditioner
*            =11: Multigrid-preconditioner
*            >50: user-defined preconditioner tag
*            This tag has to be initialized by the initialization
*            routines and to be used in the callback-routines to
*            choose a preconditioner. Identifiers 1..10 are designed
*            for standard preconditioner routines. There are predefined
*            preconditioning routines that can be called from the
*            callback routines. These are configured by parameters in the
*            PRECI/PRECD parameter block only.
*            Identifier 11 is reserved to notify the callback routine
*            that a multigrid preconditioner is used. Identifiers >50
*            are user-defined, the callback-routine must itself find
*            out how to react. In this case there is typically a
*            preconditioner-structure (i.e. a solver structure with
*            solver-tag identifying a preconditioner) or another type
*            of information array appended to this solver structure.
*            The callback-routine can then use this information it
*            correctly perform the preconditioning.
*   INRM   : Type of norm to use for convergence criteria.
*            =0: use weighted l2-norm = standard
*            =1: use euclidian vector norm
*   ITIM   : Generate timing statistics (if supported by the solver)
*            =0: don't measure time
*            =1: measure total time, put into TMTIM
*            =2: measure time of all components; timing information
*                is saved in the TMxxx-variables in the double-
*                precision structure = standard
*   IASRLN : Number of iterations that should be used to calculate
*            the asymptotic convergtence rate, if the solver supports
*            to generate an asymptotic convergence rate of the last
*            couple of steps. 0 <= IASRLN <= 32.
*
*   // Output-Tags for messages on the screen. MSGTRM is the output
*   // level for the solver itself. The other output levels are ignored
*   // normally by the algorithm itself, but should be used by the
*   // callback routines to determine their level of output.
*
*   MSGTRM : Level of output of the solver of algorithm.
*            =0: no output
*            =1: basic output
*            =2: standard output = standard
*   MFLTRM : Level of output to terminal for filter
*            =0: no output = standard
*            =1: basic output
*            =2: extended output
*   MTRMRS : Restriction of output to specific iterations.
*            Only valid for iteration-based solver algorithms.
*            Restricts the output that way, that only in each
*            MTRMRS'th iteration the residuum will be printed.
*            Standard = 1
*
*   // Standard blocks for matrix-vector multiplication and
*   // preconditioning. These blocks are reserved for standard
*   // callback-routines that need only a limited set of
*   // information.
*   // Standard MV-routines expect matrix information at position 32
*   // in the parameter array. Standard preconditioning routines
*   // expect information at position 48 in the parameter array.
*
*   MATI   : array [1..16] of integer
*            Saves general information for matrix-vector multiplication.
*            Only used, if standard matrix-vector multiplication
*            callback routines are used in the solver, which access this
*            information. Otherwise filled with 0.
*            Can be used arbitrary by user-defined matrix-vector
*            multiplication routines.
*            Not used in any black-box solver internally.
*   PRECI  : array [1..16] of integer
*            Saves general information for the preconditioner.
*            Only used, if standard preconditioner callback
*            routines are used in the solver, which access this
*            information. Otherwise filled with 0.
*            Can be used arbitrary by user-defined preconditioning
*            routines. Not used in any black-box solver internally.
*
*   // Status variables, only valid during the solution process.
*   // The solver must make sure that these values are always
*   // up-to-date!
*
*   CRITE  : Current iteration step
*   NEQ    : Number of equations in the solution/RHS vectors
*
*   // User-defineable multi-purpose block
*
*   IUDSLA : array [1..32] of double precision
*            User defined parameter block that can be customized
*            to the need of a solver. Normally filled with 0.
* end
*
* TSolverDParams = record
*
*   // Output parameters
*
*   DEFINI : initial defect
*   DEFFIN : final defect after last iteration
*   RHO    : convergence rate
*   RHOASM : asymptotic convergence rate, calculated from the last
*            IASRLN multigrid steps (if supported by the solver)
*
*   // Timing output parameters; only set if IMTIM > 0:
*
*   TMTOT  : total time for solver
*   TMFILT : time for filtering (if timing is supported by the solver)
*
*   // Input parameters
*
*   OMEGA  : General solver parameter, solver specific.
*            Can be set to 0 if not used.
*   OMGPC  : General preconditioning parameter, solver specific.
*            Can be set to 0 if not used.
*   EPSREL : relative stopping criterion. Stop iteration if
*            !!defect!! < EPSREL * !!initial defect!!
*            standard = 1E-5
*            Remark: 0=use absolute stopping criterion
*   EPSABS : absolute stopping criterion. Stop iteration if
*            !!defect!! < EPSABS
*            standard = 0 = use relative stopping criterion
*            Remark: Don't set EPSREL=EPSABS=0 !
*   DIVREL : relative divergence criterion. Treat iteration as
*            diverged if
*            !!defect!! >= DIVREL * !!initial defect!!
*            standard = 1E3
*   DIVABS : absolute divergence criterion. Treat iteration as
*            diverged if
*            !!defect!! >= DIVABS
*            standard = 1E99 = infinity
*   VECZER : RHS-vector is treated as zero im
*            max(defect) < VECZER
*            standard = 1E-12
*
*   // Standard blocks for matrix-vector multiplication and
*   // preconditioning. These blocks are reserved for standard
*   // callback-routines that need only a limited set of
*   // information.
*   // Standard MV-routines expect matrix information at position 32
*   // in the parameter array. Standard preconditioning routines
*   // expect information at position 48 in the parameter array.
*
*   MATD   : array [1..16] of integer
*            Saves general information for matrix-vector multiplication.
*            Only used, if standard matrix-vector multiplication
*            callback routines are used in the solver, that access this
*            information. Otherwise filled with 0.
*   PRECD  : array [1..16] of integer
*            Saves general information for the preconditioner.
*            Only used, if standard preconditioner callback
*            routines are used in the solver, that access this
*            information. Otherwise filled with 0.
*            This can be used to store more detailed information,
*            if a preconditioner is dependent not only on one
*            double parameter. For standard cases (only one double
*            parameter necessary), the preconditioner information
*            should be saved in OMGPC.
*
*   // User-defineable multi-purpose block
*
*   DUDSLA : array [1..32] of double precision
*            User defined parameter block that can be customized
*            to the need of a solver. Normally filled with 0.
*
* end
*
* The standard settings of the parameters can be initialised with
* INGSLV.
*
* This solver structure is generally also used for extended
* preconditioners and smoothers in a multigrid algorithm. In this
* case most variables are simply unused. Which variables of these
* blocks are used is documented at the smoother routine.
*
* EXAMPLES
* --------
* There are multiple ways to use this structure for solvers,
* preconditioners and smoothers. Here we give some examples and
* guidelines about the usage, i.e. when to use which combination
* of structure and how to set the parameters for configuring
* matrix/vector multiplication and preconditioning.
*
* 1.) Standard SSOR-solver, standard structure-7 matrix,
*     no preconditioning
*
* a) Initialize a standard solver structure
* b) Initialize the solver tag to "SSOR-solver"
* c) Initialise IMAT according to the standard matrix-vector
*    multiplication callback routines
* c) Start the SSOR-solver with this structure as parameter,
*    Standard structure-7 callback routine for MV-mult.
* d) Release the solver structure
*
* 2.) Standard CG-solver, standard structure-7 matrix,
*     standard Jacobi preconditioning
*
* a) Initialize a standard solver structure
* b) Initialize the solver tag to "CG-solver"
* c) Initialise IMAT according to the standard matrix-vector
*    multiplication callback routines
* d) Initialise PCTAG to "Jacobi preconditioning"
* e) Initialise PRECI/PRECD according to the standard Jacobi
*    preconditioning callback routines
* f) Start the CG-solver with this structure as parameter
* g) Release the solver structure
*
* 3.) Standard BiCGStab-solver, standard structure-7 matrix,
*     preconditioning dependent on variable in a DAT-file
*
* a) Initialize an extended solver structure consisting of:
*    - one standard solver block for the BiCGStab solver
*    - one solver block for the preconditioner
* b) Initialize the solver tag to "BiCGStab-solver"
* c) Initialise IMAT according to the standard matrix-vector
*    multiplication callback routines
* d) Set PCTAG=51 to indicate user-defined preconditioning
* e) Initialise the SLTAG-variable in the solver block for the
*    preconditioner to the variable in the DAT-file identifying
*    the preconditioner
* f) Write a callback-routine for the preconditioner that
*    chooses depending on OSLTAG in the preconditioner structure
*    the correct preconditioner
* f) Start the CG-solver with this structure and the user-defined
*    preconditioning routine as parameter
* g) Release the solver structure
*
* 4.) Standard BiCGStab-Solver, standard structure-7 matrix,
*     multigrid preconditioner,
*     BiCGStab-smoother in multigrid with ILU-0 preconditioner,
*     CG coarse-grid solver with SSOR-preconditioner
*
*     This is a highly advanced solver combination which requires
*     multiple blocks after each other. We give only the rough
*     guidelines here, the initialization is similar to above:
*     - standard solver-block for BiCGStab solver, PCTAG
*       initialized to "multigrid preconditioner"
*     - multigrid preconditioner block: multigrid solver-block
*       with SLTAG initialized to "multigrid preconditioner",
*       SMTAG initiaslized to "BiCGStab-smoother"
*     - BiCGStab smoother-block: BiCGStab solver-block with
*       SLTAG initialized to "BiCGStab smoother", PCTAG initialized
*       to ILU0-preconditioner, PRECI/PRECD initialized to meet the
*       parameter structure for the standard ILU0 preconditioner
*     - CG solver block for the coarse-grid solver,
*       SLTAG initialized to "CG solver", PCTAG initialized
*       to SSOR-preconditioner, PRECI/PRECD initialized to meet the
*       parameter structure for the standard SSOR preconditioner.
*     User-defined callback-routines have to be provided to call
*     the BiCGStab-smoother / CG-solver with the correct parameters.
*     The BiCGStab solver / CG-smoother expects a BiCGStab solver
*     structure - that's the reason why there are two BiCGStab
*     structures appearent: the callback-routine calls the smoother
*     with the appropriate BiCGStab-smoother structure, and the
*     CG-solver with the appropriate CG-solver structure, respectively.
************************************************************************

*****************
* TSolverIParams
*****************

      INTEGER OSTATUS,OITE  ,ONITMIN,ONITMAX,OFILT, OMATI,OPRECI,OMTRMRS
      INTEGER OPCTAG ,OINRM  ,OTIM ,OIASRLN,OMSGTRM,OMFLTRM,ONEQ,OCRITE
      INTEGER OSLTAG, OIUDSLA
      INTEGER SZSLVI

      PARAMETER (OSTATUS =  1)
      PARAMETER (OITE    =  2)
      
C 3..15: reserved
      
      PARAMETER (OSLTAG  = 16)
      PARAMETER (ONITMIN = 17)
      PARAMETER (ONITMAX = 18)
      PARAMETER (OFILT   = 19)
      PARAMETER (OPCTAG   = 20)
      PARAMETER (OINRM   = 21)
      PARAMETER (OTIM    = 22)
      PARAMETER (OIASRLN = 23)
      PARAMETER (OMSGTRM = 24)
      PARAMETER (OMFLTRM = 25)
      PARAMETER (OMTRMRS = 26)

C 26..31: reserved

      PARAMETER (OMATI   = 32)
      PARAMETER (OPRECI  = 48)

C 64..79: reserved

      PARAMETER (OCRITE  = 80)
      PARAMETER (ONEQ    = 81)
      
C 81..99: reserved

C 100..131: user-defineable, in-block
      
      PARAMETER (OIUDSLA = 100)
      
      PARAMETER (SZSLVI  = 256)
                 
C SZSLVI+1..*   : user defineable, but outside of standard-block
                 
***************
* TM020DParams
***************
                 
      INTEGER ODEFINI,ODEFFIN,ORHO   ,ORHOASM,OTMTOT ,OOMEGA,OOMGPC
      INTEGER OTMFILT,OEPSREL,OEPSABS,ODIVREL,ODIVABS,OVECZER
      INTEGER OMATD,OPRECD,ODUDSLA
      INTEGER SZSLVD
                 
      PARAMETER (ODEFINI =  1)
      PARAMETER (ODEFFIN =  2)
      PARAMETER (ORHO    =  3)
      PARAMETER (ORHOASM =  4)
      
C 5..15: reserved
      
      PARAMETER (OTMTOT  = 16)
      PARAMETER (OTMFILT = 17)

C 18..31: reserved
      
      PARAMETER (OMATD   = 32)
      PARAMETER (OPRECD  = 48)

C 64..79: reserved
      
      PARAMETER (OOMEGA  = 80)
      PARAMETER (OOMGPC  = 81)
      PARAMETER (OEPSREL = 82)
      PARAMETER (OEPSABS = 83)
      PARAMETER (ODIVREL = 84)
      PARAMETER (ODIVABS = 85)
      PARAMETER (OVECZER = 86)

C 87..99  : reserved

C 100..128  : user-defineable, in-block

      PARAMETER (ODUDSLA = 100)

      PARAMETER (SZSLVD  = 256)

C SZSLVD+1..* : user-defineable, but outside of standard-block
