************************************************************************
* This file describes the parameter structure which is passed to the
* optimization routines after being initialized by INIOPT. The
* structures are transferred to the optimizer and to the callback-
* routines of the optmizer and determine its type and behaviour.
*
* Both structures are realized by arrays
*
*   IPARAM = array [1..SZOPTI] of integer
*            = TIOptimizationParams
*   DPARAM = array [1..SZOPTD] of double
*            = TDOptimizationParams
*
* The different position in these arrays represent the parameters,
* the access can be done using the offset constants defined here,
* e.g.
*     IOPTTP = IPARAM(OIOPTTP)
*
* The content of the structures is defined as follows. Each structure
* contains three parts. The "Input"-part define pure input parameters,
* which determine the behaviour of the optimizer. The "Status"-part
* on the other hand contains dynamic information which is only
* valid during the optimization process and which is maintained
* by the optimizer and/or the callback routines. Their content before
* the beginning and after the end of the optimizer is undefined.
* The "Output" part at last contains output variables of the optimizer
* to the caller and can be used.
*
* TIOptimizationParams = record
*  
*   // Input parameters
*
*   IOPTTP - type of optimisation (manually by file / hardcoded)
*   IOPTPR - Problem identifier; can be used to identify the problem
*            which is to calculate (if the code supports this)
*   OITSKP - If > 0, the first OITSKP optimization iterations are
*            skipped. This is only used in brute-force-tests where 
*            no real optimization is performed.
*   OITMIN - Minimum number of optimization iterations; Standard=0. 
*   OITMAX - Maximum number of optimization iterations; < 0: automatic
*   IOPTCC - type of drag/lift calculation method (normal, interpolated,...)
*   IOPTFL - Type of functional to optimize for. Geometry- and 
*            problem dependent.
*   IOPTNM - Norm to use for calculation of the value of the functional. 
*   IOPTCF - Configuration flag bitfield for optimization algorithm. 
*            Algorithm-specific.
*
*   // Status variables
*
*   OPTCNT - optimization counter - stores the number of the completed
*            optimization steps
*   IXNUM  - Current X-coordinate index in brute force tests
*   IYNUM  - Current Y-coordinate index in brute force tests
*   IRNUM  - Current index of rotation coordinate in brute force tests
*
* end
*
* TDOptimizationParams = record
*  
*   // Input parameters
*
*   DVIARX - Relaxation of ALPHA-calculation in drag/lift evaluation
*   DOPEPS - Fitness parameter EPS. Optimization algorithm stops if norm of
*            functional drops below this value.
*   DOPINS - Stepsize-parameter for 1D-based optimization algorithms. 
*            Algorithm dependent. Line-search: Armijo step length control 
*            parameter DELTA, standard=0.5D0
*
*   // Input-parameters for brute-force tests
*
*   OBXMIN - minimum X-coordinate
*   OBXMAX - maximum X-coordinate
*   OBXINI - initial X-coordinate if optimization algorithm is used. 
*            For single tests, DCXPOS from geometry.dat is used instead.
*   OBXSTP - X-stepsize; Must be > 0!
*   OBYMIN - minimum Y-coordinate (0.06)
*   OBYMAX - maximum Y-coordinate (0.35)
*   OBYINI - initial Y-coordinate if optimization algorithm is used. 
*            For single tests, DCYPOS from geometry.dat is used instead.
*   OBYSTP - Y-stepsize; Must be > 0!
*   OBRMIN - minimum rotation angle 
*   OBRMAX - maximum rotation angle 
*   OBRINI - initial rotation if optimization algorithm is used. 
*            For single tests, DCROT from geometry.dat is used instead.
*   OBRSTP - rotation angle stepsize; Must be > 0!
*
*   // Status variables
*
* end
*
************************************************************************

      INTEGER OIOPTTP,OOITSKP,OOITMIN,OOITMAX,OIOPTCC,OIOPTFL,OIOPTNM
      INTEGER OIOPTCF,OOPTCNT,OIXNUM ,OIYNUM ,OIRNUM, OIOPTPR 

      PARAMETER (OIOPTTP = 1)
      PARAMETER (OIOPTPR = 2)
      PARAMETER (OOITSKP = 3)
      PARAMETER (OOITMIN = 4)
      PARAMETER (OOITMAX = 5)
      PARAMETER (OIOPTCC = 6)
      PARAMETER (OIOPTFL = 7)
      PARAMETER (OIOPTNM = 8)
      PARAMETER (OIOPTCF = 9)
      PARAMETER (OOPTCNT = 10)
      PARAMETER (OIXNUM  = 11)
      PARAMETER (OIYNUM  = 12)
      PARAMETER (OIRNUM  = 13)

      INTEGER SZOPTI
      PARAMETER (SZOPTI = 16)
      
      INTEGER ODVIARX,ODOPEPS,ODOPINS,OOBXMIN
      INTEGER OOBXMAX,OOBXINI,OOBXSTP,OOBYMIN,OOBYMAX,OOBYINI,OOBYSTP
      INTEGER OOBRMIN,OOBRMAX,OOBRINI,OOBRSTP

      PARAMETER (ODVIARX = 1)
      PARAMETER (ODOPEPS = 2)
      PARAMETER (ODOPINS = 3)
      PARAMETER (OOBXMIN = 4)
      PARAMETER (OOBXMAX = 5)
      PARAMETER (OOBXINI = 6)
      PARAMETER (OOBXSTP = 7)
      PARAMETER (OOBYMIN = 8)
      PARAMETER (OOBYMAX = 9)
      PARAMETER (OOBYINI = 10)
      PARAMETER (OOBYSTP = 11)
      PARAMETER (OOBRMIN = 12)
      PARAMETER (OOBRMAX = 13)
      PARAMETER (OOBRINI = 14)
      PARAMETER (OOBRSTP = 15)

      INTEGER SZOPTD
      PARAMETER (SZOPTD = 16)