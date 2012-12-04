************************************************************************
* This files describes the solver parameter structure of the grid
* adaption routines. This parameter structure contains input
* parameters that configure how the linear solver in the grid
* adaption behaves. There are no output variables in these
* structures, nor are there any status variables.
*
* The realization of the structures as arrays is made as:
*
*   IPARAM = array [1..SZGLSI] of integer 
*          = TIGridAdaptionSolver
*   DPARAM = array [1..SZGLSD] of double
*          = TDGridAdaptionSolver
*
* The Oxxxxx-constants below can be used to access the different
* variables inside of the structure-array, e.g.:
*
*   IGSSLV = IPARAM(OIGSSLV)
*
* TIGridAdaptionSolver = record
*   // Input parameters:
*   IGSSLV : Type of solver for grid adaption
*            =0: BiCGStab
*            =1: Multigrid
*            =2: UMFPACK4
*   ICUBGS : Cubature formula for building grid adaption matrix
*   ICURGS : Cubature formula for building grid adaption RHS
*   IRCGRD : Type of method to use for reconstructing gradients.
*            =0: ZZ-technique (higher order)
*            =1: Standard interpolation technique (1st order)
*   IRELGS : Relative or absolute error bound for solving system
*            =0: absolute
*            =1: relative
*   NITSGS : Maximum number of solver iterations for grid adaption solver
*   IPCSMG : Type of smoother/preconditioner; 
*            =-1: Jac/SSOR (depending on OMGPGS, see II010)
*            =0=1=2=3: BiCGStab(IPCSGS)
*   NSMSGS : MG-Solver: Number of pre-/postsmoothing steps if MG is 
*            used as solver
*   IPCSGS : MG-Solver: Type of preconditioner if BiCGStab-Smoother is chosen
*            =0: Jacobi
*            =1: SSOR(OMGPGS)
*            =2: ILU(0)
*   ISRTGS : MG-Solver: Type of resorting
*            =0: none
*            =1: Cuthill-McKee
*   ICGSGS : MG-Solver: Type of coarse grid solver
*            =0: BiCGStab
*            =1: UMFPACK4
*   ICGSLV : MG-Solver: Minimum level, level where the coarse-grid 
*            solver solves
*   NICGGS : MG-Solver: Maximum number of iterations in coarse 
*            grid solver
*   IRCGGS : MG-Solver: Use relative or absolute error bound for coarse
*            grid solver
*            =0: absolute
*            =1: relative
*   IPCGGS : MG-Solver: Type of preconditioner for coarse grid solver
*            =-1: Jac/SSOR (depending on OMCGGS, see II010)
*            =0: Jacobi
*            =1: SSOR(OMGPGS)
*            =2: ILU(0)
* end
* 
* TDGridAdaptionSolver = record
*   // Input parameters:
*   EPSGSS : Accuracy of grid-apation solver
*   OMGPGS : Damping parameter of smoother/preconditioner
*   EPSCGS : MG-Solver: Accuracy of coarse grid solver
*   OMCGGS : MG-Solver: Damping parameter of preconditioner in coarse
*            grid solver
* end
************************************************************************

***********************
* TIGridAdaptionSolver
***********************

      INTEGER OIGSSLV,OICUBGS,OICURGS,OIRELGS,ONITSGS,OIPCSMG,ONSMSGS
      INTEGER OIPCSGS,OISRTGS,OICGSGS,OICGSLV,ONICGGS,OIRCGGS,OIPCGGS
      INTEGER OIRCGRD

      PARAMETER (OIGSSLV=1)
      PARAMETER (OICUBGS=2)
      PARAMETER (OICURGS=3)
      PARAMETER (OIRCGRD=4)
      PARAMETER (OIRELGS=5)
      PARAMETER (ONITSGS=6)
      PARAMETER (OIPCSMG=7)
      PARAMETER (ONSMSGS=8)
      PARAMETER (OIPCSGS=9)
      PARAMETER (OISRTGS=10)
      PARAMETER (OICGSGS=11)
      PARAMETER (OICGSLV=12)
      PARAMETER (ONICGGS=13)
      PARAMETER (OIRCGGS=14)
      PARAMETER (OIPCGGS=15)

C 15..24: reserved
C 25..31: user-defineable

      INTEGER SZGLSI
      PARAMETER (SZGLSI=32)
      
***********************
* TDGridAdaptionSolver
***********************

      INTEGER OEPSGSS,OOMGPGS,OEPSCGS,OOMCGGS
      
      PARAMETER (OEPSGSS=1)
      PARAMETER (OOMGPGS=2)
      PARAMETER (OEPSCGS=3)
      PARAMETER (OOMCGGS=4)

C 5..24: reserved
C 25..31: user-defineable

      INTEGER SZGLSD
      PARAMETER (SZGLSD=32)
