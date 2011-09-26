************************************************************************
* This file describes the solver structure for the nonlinear
* iteration for solving the stationary Navier Stokes equation (NSDEF2).
* As the structure is an extension to the general solver structure,
* the file "SSOLVERS.INC" has to be included before.
*
* It defines two parameter blocks for integer and double parameters,
* respectively. Both parameter blocks contain three parts:
* - An input-block which has to be initialized by the caller
* - An output-block which returns the result of the computation
* - A status- & temporary block, which is only valid
*   during the computation.
* The caller can attach additional elements to the parameter blocks
* with further information for the callback-routines, e.g. matrix-
* handles, flags, tags,...
*
* The structures themselves are realized as integer/double precision
* arrays:
*   IPARAM = array [1..SZNSDI] of integer 
*          = TNSDefIParams
*   DPARAM = array [1..SZNSDD] of double
*          = TNSDefDParams
*
* These structures are direct extensions to the original solver
* structure, i.e. there has no SZxxxx-constant to be added to
* the offset variable, the variables can be acessed directly
* by a construction like IASMBL(OIALPHA).
* The structures are organized as follows. The first part of the
* structure is as usual a block with general solver parameters.
* Some of the general parameters are not used as indicated with "=0"
* below.
* Attached to that there is a couple of solver specific variables
* for NSDEF as well as a definition block for the enclosed
* multigrid solver for the linear sub-problems. Therefore, the
* include files SLINSOL.INC, SLINSOL2DNS.INC and STRIA.INC have to be 
* included prior to this file.
* The structure contains timing information, therefore the include
* file STIMING.INC has to be included before.
* At the end of the structure, a substructure for the multigrid
* solver is embedded. Therefore, the file M020.INC has to be
* included before.
*
* The whole structure can be initializes with the INISTS subroutine
* in NSDEF.F.
*
* TNSDIParams = record 
*   // The first part of TNSDIParams is the whole TSolverIParams-
*   // structure, i.e. TNSDIParams is a direct extension of
*   // TSolverIParams:
*
*   // Standard-Output:
*
*   STATUS : Result of the solving process.
*            =0: Solver converged successfully
*            =1: Nonlinear solver broke down
*            =2: Linear solver broke down
*   ITE    : Number of performed iterations
*
*   // Standard-Input:
*
*   SLTAG  : Solver tag.
*            =20: Stationary Navier-Stokes, coupled method.
*   NITMIN : minimum number of iterations
*   NITMAX : maximum number of iterations
*   IFILT  : = 0
*   OPCTAG : = 0
*   INRM   : = 0
*   ITIM   : = 0
*   IASRLN : = 0
*   MSGTRM : Level of output of the solver of algorithm; standard = 2
*   MFLTRM : = 0
*   MTRMRS : = 1
*   MATI   : = 0
*   PRECI  : = 0
*
*   // Standard-Status
*
*   CRITE  : Current iteration step
*   NEQ    : Number of equations in the solution/RHS vectors
*
*   // The second part now describes the configuration of the
*   // nonlinear solver specific settings. It again consists of an
*   // Output-block, an Input-block and a Status-block.
*
*   // Nonlinear solver, Output:
*
*   NLIN    : Total number of iterations of the linear solver
*
*   // Nonlinear solver, Input:
*
*   MSGLIN : Message level for solver of linear equations (multigrid)
*   
*   IAUSAV : Whether or not to save the intermediate solution of the
*            nonlinear iteration to disc.
*            = 0: Don't save the solution.
*            > 0: Save the solution every IAUSAV steps to a file
*                 data/AUTOSAV.xxxx
*
*   LFLAUT : Handle to a DString with the basic path/filename
*            for saving intermediate solutions to disc if IAUSAV<>0.
*
*   // Nonlinear solver, Status:
*
* end
*
* TNSDDParams = record 
*   // The first part of TNSDDParams is the whole TSolverDParams-
*   // structure, i.e. TNSDDParams is a direct extension of
*   // TSolverDParams:
*
*   // Standard-Output:
*
*   DEFINI : = 0
*   DEFFIN : = 0
*   RHO    : = 0
*   RHOASM : = 0
*
*   TMTOT  : total time for nonlinear solver
*   TMFILT : = 0
*
*   // Standard-Input:
*
*   OMEGA  : = 0
*   OMGPC  : = 0
*
*   // The standard stopping criteria of the solver structure are not
*   // used! Instead, the structure defines a couple of "specialized"
*   // stopping criteria for the nonlinear solver, which cannot be
*   // acquires by the standard stopping criteria - see below.
*
*   EPSREL : = 0
*   EPSABS : = 0
*   DIVREL : = 0
*   DIVABS : = 0
*
*   VECZER : = 0
*   MATD   : = 0
*   PRECD  : = 0
*
*   // Standard-Status
*
*
*   // The second part now describes the configuration of the
*   // nonlinear solver specific settings. It again consists of an
*   // Output-block, an Input-block and a Status-block.
*
*   // Nonlinear solver, Output:
*
*   // Timing information
*
*   TNLTIM : array [1..SZTIMG] of double
*            TTimings structure: Timing information of the
*            different components of the algorithm.
*
*   // Nonlinear solver, Input:
*
*   OMGINI : initial OMEGA parameter for optimal defect correction
*   OMGMIN : minimum OMEGA parameter for optimal defect correction
*   OMGMAX : maximum OMEGA parameter for optimal defect correction
*
*   // The next block of parameters defines a couple of "specialized" 
*   // stopping criteria for the nonlinear iteration. They replace the
*   // standard stopping criteria above!
*
*   EPSD   : limit for U-defects
*   EPSDIV : limit for divergence-defects
*   EPSUR  : limit for relative U-changes
*   EPSPR  : limit for relative P-changes
*   DMPD   : limit for defect improvement
* 
*   // Nonlinear solver, Status:
*
* end
************************************************************************

      INTEGER ONLIN  ,OMSGLIN
      INTEGER OIAUSAV,OLFLAUT

      INTEGER SZNSDI
      
      PARAMETER (ONLIN   = SZSLVI +  1)
      PARAMETER (OMSGLIN = SZSLVI +  2)
      PARAMETER (OIAUSAV = SZSLVI +  3)
      PARAMETER (OLFLAUT = SZSLVI +  4)

      PARAMETER (SZNSDI  = OIAUSAV + 32)
      
      INTEGER SZNSDD
      
      
      INTEGER OTNLTIM
      INTEGER OOMGINI,OOMGMIN,OOMGMAX
      
      INTEGER OEPSD  ,OEPSDIV,OEPSUR ,OEPSPR ,ODMPD 
      
      PARAMETER (OTNLTIM = SZSLVD +  1)
      PARAMETER (OOMGINI = SZSLVD + SZTIMG +  1)
      PARAMETER (OEPSD   = SZSLVD + SZTIMG +  2)
      PARAMETER (OEPSDIV = SZSLVD + SZTIMG +  3)
      PARAMETER (OEPSUR  = SZSLVD + SZTIMG +  4)
      PARAMETER (OEPSPR  = SZSLVD + SZTIMG +  5)
      PARAMETER (ODMPD   = SZSLVD + SZTIMG +  6)
      PARAMETER (OOMGMIN = SZSLVD + SZTIMG +  7)
      PARAMETER (OOMGMAX = SZSLVD + SZTIMG +  8)
      
      PARAMETER (SZNSDD  = OOMGMAX + 32)
      