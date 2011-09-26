************************************************************************
* This file describes the internal variable representation block
* used in the solver of the static grid adaption. This block is 
* passed to the solver, and from the solver to the different
* user-defined subroutines. It contains the information about the
* matrix handles,... which allowes all routines to work without
* the use of COMMON blocks.
*
* In contrast to the dynamic grid adaption routines, the static grid
* adaption routines only use the extended multigrid solver M020
* and reduce the usage of COMMON-blocks. The basis for
* solving the arising system is a pair of TM020IParams/TM020DParams
* structure, extended by a couple of variables for storing
* information about the iteration, and a set of callback routines
* for the M020. 
*
* There are two structures defined here: TM020IExtParams and
* TM020DExtParams for integer and double-precision variables.
* Both structures enhance the original M020 multigrid solver
* structure, i.e. they are direct extensions to it. This means
* both structures contain the original TM020IParams and
* TM020DParams structures at the beginning and add additional
* variables to the end. 
*
* The access to the content can be done using the offset constants in 
* this file, e.g.:
*
*   ISRT = IPARAM(OISRT)
*
* These structures are direct extensions to the original solver
* structure, i.e. there has no SZxxxx-constant to be added to
* the offset variable.
*
* The structures are realized as arrays:
*   IPARAM = array [1..SZGSSI] of integer 
*          = TIGridAdaptionCommon
*   DPARAM = array [1..SZGSSD] of double
*          = TDGridAdaptionCommon
* The mapping of the structure to the array can generally be seen 
* as follows:
*
* TM020IExtParams 
*   Array-Position   Function
*   =========================
*    1         --
*    ...        |--- General solver structure, input/output variables
*    SZSLVI    --    defining the general behaviour; see SSOLVERS.INC
*
*    SZSLVI+1  --    
*    ...        |--- Extended solver structure, input/output variables
*    SZ020I    --    defining the behaviour of the M020 solver; M020.INC
*
*    SZ020I+1  --    Specializes solver structure with information for
*    ...        |--- the callback-routines of M020; this is defined 
*    SZGSSI    --    here; handles and parameters; only input-variables
*                    
*
* The same structure holds for TM020DExtParams.
*
* For a proper definition of the constants here, the include files 
* SSOLVERS.INC and M020.inc have to be included prior to this file, 
* as well as the file CBASICMG.INC, which is the basis for M020.INC,
* and SZTRIA, as triangulation structures are saved in the multigrid
* structure.
* The detailed definition of the structures is as follows.
* Data that is not used in a special configuration is set to 0.
*
* TM020IExtParams = record
*
*   // The first part of TM020IParams is the whole TM020IParams-
*   // structure, which itself include the TSolverIParams solver
*   // structure:
*
*   STATUS : Result of the solving process.
*   ITE    : Number of performed multigrid iterations
*   ...
*   NEQ    : Number of equations in the solution/RHS vectors
*   ...
*   KIT0   : array [1..NNLEV] of integer
*
*   // Additional input-parameters, have to be initialized
*   // before calling the solver and released when the solver
*   // is finished. Are used in the callback-routines
*   // for correct operation.
*  
*   // Handle-Block: Contains handles of matrices for the system.
*
*   LA      : array [1..NNLEV] of integer
*             Handles of the matrix entries of system matrices on all levels
*   LCOL    : array [1..NNLEV] of integer
*             Handles of the matrix column structure on all levels
*   LLD     : array [1..NNLEV] of integer
*             Handles of the matrix row structure on all levels
*   LSOL    : array [1..NNLEV] of integer
*             Handles of the solution vectors on all levels
*   LRHS    : array [1..NNLEV] of integer
*             Handles of the RHS vectors on all levels
*   LAUXD   : array [1..NNLEV] of integer
*             Handles of auxiliary vectors on all levels
*   LAUX2   : array [1..NNLEV] of integer
*             Handles of auxiliary vectors on all levels
*   KNA     : array [1..NNLEV] of integer
*             Size of the matrices on all levels
*   KNEQ    : array [1..NNLEV] of integer
*             Number of equations on all levels
*   ISRT    : integer
*             Determines the currently active sorting strategy.
*             = 0: matrices are unsorted, no sorting strategy has been set up
*             =-1: matrices are unsorted, Cuthill McKee sorting
*                  strategy has been set up.
*             = 1: matrices are sorted for Cuthill-McKee,
*                  permutation and inverse permutation are to be found
*                  in KTR1/KTR2
*   KTR1    : array [1..NNLEV] of integer
*             For each level, handle to a permutation, how the vectors
*             must be permuted to fit to the matrix
*   KTR2    : array [1..NNLEV] of integer
*             Handle to the inverse permutation to KTR1(.) on
*             each level
*   KILU0   : array [1..NNLEV] of integer
*             Handle to entries of ILU0-matrices on all levels
*   KU4SYM  : Handle to UMFPACK4 symbolical factorization on 
*             coarse grid
*   KU4NUM  : Handle to UMFPACK4 numerical factorization on coarse grid
*   SMSI    : TSolverIParams
*             Structure to configure the smoother
*   CGSI    : TSolverIParams
*             Structure to configure the coarse grid solver
*
*   // Triangulation block: contains the triangulation structures
*   // of all levels. Necessary for prolongation/restruction
* 
*   TRIAS   : array [1..NNLEV] of STRIA
*             = array [1..SZTRIA,1..NNLEV] of integer
*             Triangulation structures on all levels
* end
*
* TM020DExtParams = record
*
*   // The first part of TM020DParams is the whole TM020DParams-
*   // structure, which itself include the TSolverDParams solver
*   // structure:
*
*   DEFINI : initial defect
*   DEFFIN : final defect after last iteration
*   ...
*   PRECD  : array [1..16] of integer
*   ...
*   TMCGC  : time for coarse-grid corrections
*   
*   // Additional input-parameters, have to be initialized
*   // before calling the solver and released when the solver
*   // is finished.  Are used in the callback-routines
*   // for correct operation.
*
*   CNTRU4  : array [1..20] of double
*             control structure for UMFPACK4 coarse grid solver,
*             if used.
*   INFOU4  : array [1..90] of double
*             info structure for UMFPACK4 coarse grid solver,
*             if used
*   SMSD    : TSolverDParams
*             Structure to configure the smoother.
*             This contains a complete solver structure which
*             can be prepared for the use as a smoother
*             if necessary. The output variables here are only used
*             in the callback-routine.
*   CGSD    : TSolverDParams
*             Structure to configure the coarse grid solver.
*             This contains a complete solver structure which
*             can be prepared for the use as a coarse-grid solver
*             if necessary. The output variables here are only used
*             in the callback-routine.
* end
*
************************************************************************

******************                                      
* TM020IExtParams
******************

      INTEGER OLA    ,OLCOL  ,OLLD   ,OLSOL  ,OLRHS  ,OLAUXD,OLAUX2,OKNA   
      INTEGER OISRT  ,OLTR1  ,OLTR2  ,OLILU0 ,OKU4SYM,OKU4NUM
      INTEGER OCGSI, OSMSI, OTRIAS

      INTEGER SZGSSI,SZGSSD

      PARAMETER (OLA     = SZ020I+1      )
      PARAMETER (OLCOL   = OLA    +NNLEV )
      PARAMETER (OLLD    = OLCOL  +NNLEV )
      PARAMETER (OLSOL   = OLLD   +NNLEV )
      PARAMETER (OLRHS   = OLSOL  +NNLEV )
      PARAMETER (OLAUXD  = OLRHS  +NNLEV )
      PARAMETER (OLAUX2  = OLAUXD +NNLEV )
      PARAMETER (OKNA    = OLAUX2 +NNLEV )
      PARAMETER (OISRT   = OKNA   +NNLEV )
      PARAMETER (OLTR1   = OISRT  +1     )
      PARAMETER (OLTR2   = OLTR1  +NNLEV )
      PARAMETER (OLILU0  = OLTR2  +NNLEV )
      PARAMETER (OKU4SYM = OLILU0 +NNLEV )
      PARAMETER (OKU4NUM = OKU4SYM+1     )
      PARAMETER (OSMSI   = OKU4NUM+1     )
      PARAMETER (OCGSI   = OSMSI + SZSLVI)
      PARAMETER (OTRIAS  = OCGSI + SZSLVI)

C OCGSI + SZSLVI .. SZGSSI: reserved

      PARAMETER (SZGSSI  = OTRIAS + NNLEV*SZTRIA) 
  
******************                                      
* TM020DExtParams
******************

      INTEGER OCNTRU4,OINFOU4
      INTEGER OCGSD, OSMSD

      PARAMETER (OCNTRU4 = SZ020D +1     )
      PARAMETER (OINFOU4 = OCNTRU4 +20   )
      PARAMETER (OSMSD   = OINFOU4 +90   )
      PARAMETER (OCGSD   = OSMSD +SZSLVD )

C OCGSD + SZSLVD .. SZGSSD: reserved

      PARAMETER (SZGSSD  = OCGSD + SZSLVD +16)
