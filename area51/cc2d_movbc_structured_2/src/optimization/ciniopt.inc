************************************************************************
* This file contains the variable representing the user provided
* information in the DAT file of the application.
* More specifically the variables here define the behavior of
* the optimiser, if active.
* The data is filled by the routine RDOPT by the data of the
* optimization DAT file. Using the routine INIOPT, the data is
* transferred into a local optimization structure (c.f. 
* soptimization.inc), which is passed to the optimizer and which
* determines its behaviour.
************************************************************************

      INTEGER IOPTTP,OPTCNT,OITSKP,OITMIN,OITMAX,IOPTCC
      INTEGER IXNUM, IYNUM, IRNUM
      INTEGER IOPTFL, IOPTNM, IOPTCF, IOPTPR
      COMMON /OPTINT/ IOPTTP,OPTCNT,OITSKP,OITMIN,OITMAX,
     *                IOPTCC,IOPTPR,
     *                IXNUM, IYNUM, IRNUM, IOPTFL, IOPTNM, IOPTCF
      
      DOUBLE PRECISION DVIARX,DOPEPS,DOPINS,GEOTIM
      COMMON /OPTDBL/ DVIARX,DOPEPS,DOPINS,GEOTIM
           
      SAVE /OPTINT/, /OPTDBL/
      
C IOPTTP - type of optimisation (manually by file / hardcoded)
C IOPTPR - Identifier to choose the problem to optimize
C OPTCNT - optimization counter - stores the number of the completed
C          optimization steps
C OITSKP - If > 0, the first OITSKP optimization iterations are
C          skipped. This is only used in brute-force-tests where 
C          no real optimization is performed.
C OITMIN - Minimum number of optimization iterations; Standard=0. 
C OITMAX - Maximum number of optimization iterations; < 0: automatic
C IOPTCC - type of drag/lift calculation method (normal, interpolated,...)
C DVIARX - Relaxation of ALPHA-calculation in drag/lift evaluation

C IOPTFL  - Type of functional to optimize for. Geometry- and 
C           problem dependent.
C IOPTNM  - Norm to use for calculation of the value of the functional. 
C IOPTCF  - Configuration flag bitfield for optimization algorithm. 
C           Algorithm-specific.
C DOPEPS  - Fitness parameter EPS. Optimization algorithm stops if norm of 
C           functional drops below this value.

      DOUBLE PRECISION OBXMIN,OBXMAX,OBXSTP,OBYMIN,OBYMAX,OBYSTP,
     *                 OBRMIN,OBRMAX,OBRSTP,OBXINI,OBYINI,OBRINI
      COMMON /OPTBRT/  OBXMIN,OBXMAX,OBXSTP,OBYMIN,OBYMAX,OBYSTP,
     *                 OBRMIN,OBRMAX,OBRSTP,OBXINI,OBYINI,OBRINI
      
C Variables for brute-force tests:
C
C OBXMIN - minimum X-coordinate for brute-force tests
C OBXMAX - maximum X-coordinate for brute-force tests
C OBXSTP - X-stepsize for brute-force tests
C OBYMIN - minimum Y-coordinate for brute-force tests
C OBYMAX - maximum Y-coordinate for brute-force tests
C OBYSTP - Y-stepsize for brute-force tests
C OBXINI - Initial X-position
C OBYINI - Initial Y-position
C OBRINI - Initial rotation

C Status variables; these values save the current status of the
C optimiser:
C
C IXNUM   - Current X-coordinate index in brute force tests
C IYNUM   - Current Y-coordinate index in brute force tests
C IRNUM   - Current index of rotation coordinate in brute force tests
C -> These indices are used to build the filename for GMV.
C    They are incremented for each new configuration
