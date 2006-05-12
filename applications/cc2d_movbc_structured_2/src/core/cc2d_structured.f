************************************************************************
* CC2D_MOVBC_STRUCTURED
*
* Experimental code, based on CC2D_MOVBC.
* Solves the stationary and nonstationary (low-RE-number) Navier-
* Stokes equation.
*
* COMMON-block-less implementation in all solver components - only
* the initialisation routines and the core-routines know about
* COMMON blocks, as they are only used to store the variables of the DAT
* files.
*
* Read the file INIT2.F for implementational details about how to start
* the solvers. Read the different .INC include files for a documentation
* about the structures used.
************************************************************************

      PROGRAM CC2D_MOVBC  

      USE testcode

      IMPLICIT NONE

      INCLUDE 'cmem.inc'
      
      INCLUDE 'ctimediscr.inc'
      
      INCLUDE 'sinigeometry.inc'
      
      INCLUDE 'cout.inc'
      INCLUDE 'cfileout.inc'
      
      INTEGER MSHOW
      
      INTEGER IGEOM(SZGEOI)
      DOUBLE PRECISION DGEOM(SZGEOD)
      
      INTEGER I,S
      CHARACTER CSTR*(255)
      
C     At the very beginning initialise the pseudodynamic memory 
C     management. Nothing must be done before that!!!

      CALL ZINIT(NNWORK,'feat.msg','data/cc2d.err','data/cc2d.prt',
     *                             'data/cc2d.sys','data/cc2d.trc') 

      CALL MKTEST
      STOP

C     Initialize the output channels - also at the very beginning
C     before anything is done!

      CALL OINIT (MSHOW)

C     Read DAT files, initialize general data, create parametrization.
C     this initialises IGEOM and DGEOM.

      CALL INMAIN (MSHOW,IGEOM,DGEOM)

C     Call the stationary or instationary solver

      IF (ISTAT.EQ.0) THEN
        CALL CC2DST (MSHOW,IGEOM,DGEOM)
      ELSE
        CALL CC2DNS (MSHOW,IGEOM,DGEOM)
      END IF
      
C     General cleanup

      CALL DNMAIN (IGEOM,DGEOM)

C     Print checksum of L()-array. Print memory usage. 
C     If these values are <> 0, there is definitely
C     a memory leak in the program!!!

      S = 0
      DO I=1,NNARR
        S = S+L(I)
      END DO
      IF ((S.NE.0).OR.(IWORK.NE.0)) THEN
        WRITE (CSTR,'(A)') 'Warning: Possible memory leak!'
        CALL CNOUTS (MT,MTERM,MFILE,.TRUE.,CSTR)
        WRITE (CSTR,'(A,I10)') 'L()-checksum (should be = 0) = ',S
        CALL CNOUTS (MT,MTERM,MFILE,.TRUE.,CSTR)
        WRITE (CSTR,'(A,I10)') 'Final IWORK  (should be = 0) = ',IWORK
        CALL CNOUTS (MT,MTERM,MFILE,.TRUE.,CSTR)
      END IF

1     FORMAT(79('-'))
      
      END       

************************************************************************
* CC2D - start stationary solver for standard problem
*
* In:
*   MSHOW  - Output level during initialisation phase
*   IGEOM  - Integer parameter block with geometry information
*   DGEOM  - Double precision parameter block with geometry information
*
* Out:
*   No return parameters, all output written to files/terminal.
************************************************************************

      SUBROUTINE CC2DST (MSHOW,IGEOM,DGEOM)
      
      IMPLICIT NONE

      INCLUDE 'cmem.inc'
      INCLUDE 'cout.inc'

      INCLUDE 'stria.inc'
      INCLUDE 'cbasicmg.inc'
      
      INCLUDE 'slinsol.inc'
      INCLUDE 'slinsol2dns.inc'
      INCLUDE 'smat2dns.inc'
      INCLUDE 'ssolvers.inc'
      INCLUDE 'stiming.inc'
      INCLUDE 'm020.inc'
      INCLUDE 'snsdef.inc'
      INCLUDE 'sassembly.inc'
      
      INCLUDE 'stracking.inc'
      
      INCLUDE 'sgridadaptmgstaticparams.inc'
      INCLUDE 'cgridadapt.inc'
      
      INCLUDE 'dstrings.inc'
      
      INCLUDE 'cparametrization.inc'
      INCLUDE 'cdiscretization.inc'
      INCLUDE 'cpostproc.inc'
      INCLUDE 'clinsol_cc2d.inc'

      INTEGER MSHOW
      
      INTEGER TRIAS(SZTRIA,NNLEV)
      INTEGER MATDAT(SZN2MI,NNLEV)
      INTEGER VECDAT(SZN2VI,NNLEV)

      INTEGER ISTPAR(SZNSDI)
      DOUBLE PRECISION DSTPAR(SZNSDD)
      
C     Multigrid structures provide space for one Solver, 
C     MG-solver/preconditioner, Smoother and coarse grid solver:
      
      INTEGER IMGPAR (SZ020I+3*SZSLVI)
      DOUBLE PRECISION DMGPAR (SZ020D+3*SZSLVD)
      
      INTEGER IASMBL(SZASMI)
      DOUBLE PRECISION DASMBL(SZASMD)
      
      INTEGER IGEOM(*)
      DOUBLE PRECISION DGEOM(*)
      
      INTEGER ITRACK(SZTRKI),LFNAMS(6)
      
      INTEGER LUP,LRHS,IBDTYP,LAUX
      
      EXTERNAL GNRHSV
      
      EXTERNAL I000,DFPSTD
      
C     local variables

      INTEGER I,NUVP
      DOUBLE PRECISION WCTIM,WCTINI,WCTTOT
      
C     Store the wall clock time into WCTIM. Is used to determine the
C     total time at the end of the algorithm. 

      WCTTOT = 0D0
      CALL GTMAUX (WCTIM,WCTTOT,1,0)

C     Initialize everything for the solution process:
C
C     Initialize triangulations.
C     Correct NLMIN/NLMAX into bounds 1..9 if they are out of bounds.
      
      CALL INMTRI (MSHOW,TRIAS,NLMIN,NLMAX,0,CMESH)
      
C     Initialize stationary solver,
C     Allocate memory for matrices, vectors,...

      CALL INSTSL (NLMIN,NLMAX,TRIAS,
     *             MATDAT,VECDAT,ISTPAR,DSTPAR,IMGPAR,DMGPAR,
     *             IASMBL,DASMBL,LUP,LRHS,NUVP)
     
C     Adapt the grid if activated in the user

      IF (IGAMET.EQ.3) THEN
        CALL XSMMGW (TRIAS,NLMIN,NLMAX,IGEOM,DGEOM,0,0D0)
      END IF
      
C     Initialize/Generate matrices
     
      CALL GNMATD (NLMIN,NLMAX,TRIAS,
     *             MATDAT,ISTPAR,DSTPAR,IASMBL,DASMBL)
     
C     Generate boundary/geometry information for the discretization.
C     As we are in the stationary case, we implement boundary 
C     conditions directly into the triangulation.
     
      CALL GNBDGE (NLMIN,NLMAX,TRIAS,
     *             IASMBL,DASMBL,IGEOM,DGEOM)

      IBDTYP = 0
      IF (IBDR.EQ.1) IBDTYP = 1
      IBDTYP = IOR(4+2,IBDTYP)
     
C     Generate RHS/solution vector data
     
      CALL GNSTVC (NLMIN,NLMAX,TRIAS,
     *             VECDAT,ISTPAR,DSTPAR,IASMBL,DASMBL,
     *             IGEOM,DGEOM,
     *             IBDTYP,LUP,LRHS)

C     Read the starting vector if necessary

      IF (ISTART.NE.0) THEN
        CALL RDSTRT (NLMIN,NLMAX,VECDAT,TRIAS,DWORK(L(LUP)),
     *               IINT,IAPR,IAVPR,DPREP,
     *               ISTART,CSTART)
      END IF
      
C     Preparation completed; stop the time necessary for the
C     preprocessing:

      CALL GTMAUX (WCTIM,WCTINI,1,1)

C     Call the solver

      CALL NSDEF2 (NLMIN,NLMAX,TRIAS,
     *             MATDAT,VECDAT,ISTPAR,DSTPAR,
     *             IMGPAR,DMGPAR,IASMBL,DASMBL,
     *             NUVP,
     *             DWORK(L(LUP)),DWORK(L(LRHS)))

C     Postprocessing:
C
C     Open tracking files
      
      CALL INITRK (ITRACK,-1)
      
C     Call the postprocessing;

      WRITE (MTERM,'(A)') 'Starting postprocessing...'

C     Reserve memory for an auxiliary vector beforehand.

      CALL ZNEW (NUVP,-1,LAUX,'DAUX  ')
      
C     Prepare file output
      
      I=0
      IF (ISOL.EQ.1) I=-1
      IF (ISOL.EQ.2) I=1
      
C     Fetch names of output files into DStrings

      LFNAMS (1) = STNEWC (.TRUE.,CSOL)
      LFNAMS (2) = STNEWC (.TRUE.,CFLDU )
      LFNAMS (3) = STNEWC (.TRUE.,CFLDV )
      LFNAMS (4) = STNEWC (.TRUE.,CFLDP )
      LFNAMS (5) = STNEWC (.TRUE.,CFLISO)
      LFNAMS (6) = STNEWC (.TRUE.,CFLGMV)
      
      CALL DFPSTA (NLMIN,NLMAX,
     *             TRIAS,MATDAT,VECDAT,
     *             ISTPAR,DSTPAR,
     *             IMGPAR,DMGPAR,IASMBL,DASMBL,IGEOM,DGEOM,
     *             ITRACK,NUVP,DWORK(L(LUP)),
     *             DWORK(L(LRHS)),DWORK(L(LAUX)),
     *             IGMV,IERANA,1,4,1,
     *             I*NLMAX,IFUSAV,IFPSAV,IFXSAV,LFNAMS)

C     Release memory for postprocessing
      
      DO I=6,1,-1
        CALL STDIS(LFNAMS (I))
      END DO

      CALL ZDISP (0,LAUX,'DAUX  ')
      
C     Close tracking files

      CALL DONTRK (ITRACK)     

C     Release memory, finish

      CALL DNBDGE (NLMIN,NLMAX,
     *             TRIAS,IASMBL,DASMBL)
      CALL DNSTSL (NLMIN,NLMAX,
     *             MATDAT,VECDAT,ISTPAR,DSTPAR,LUP,LRHS)
      CALL DNMTRI (NLMIN,NLMAX,TRIAS)
      
C     Stop the complete time including the initialisation
C     and postprocessing in WCTTOT.

      CALL GTMAUX (WCTIM,WCTTOT,1,1)

C     Print statistical data / timing results

      CALL STATST (MSHOW,WCTINI,WCTTOT,ISTPAR,DSTPAR)
      
1     FORMAT(79('-'))

      END
      
************************************************************************
* CC2D - start nonstationary solver for standard problem
*
* In:
*   MSHOW  - Output level during initialisation phase
*   IGEOM  - Integer parameter block with geometry information
*   DGEOM  - Double precision parameter block with geometry information
*
* Out:
*   No return parameters, all output written to files/terminal.
************************************************************************
      
      SUBROUTINE CC2DNS (MSHOW,IGEOM,DGEOM)
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'

      INCLUDE 'stria.inc'
      INCLUDE 'cbasicmg.inc'
      
      INCLUDE 'slinsol.inc'
      INCLUDE 'slinsol2dns.inc'
      INCLUDE 'smat2dns.inc'
      INCLUDE 'ssolvers.inc'
      INCLUDE 'stiming.inc'
      INCLUDE 'm020.inc'
      INCLUDE 'snsdef.inc'
      INCLUDE 'snonstead.inc'
      INCLUDE 'sassembly.inc'
      
      INCLUDE 'sinigeometry.inc'
      
      INCLUDE 'sadtstep.inc'
      
      INCLUDE 'stracking.inc'
      
      INCLUDE 'sgridadaptmgstaticparams.inc'
      INCLUDE 'cgridadapt.inc'
      
      INCLUDE 'cdiscretization.inc'
      INCLUDE 'cparametrization.inc'
      INCLUDE 'clinsol_cc2d.inc'
      
      INCLUDE 'snonsteadpostproc.inc'

      INTEGER MSHOW
      
      INTEGER TRIAS(SZTRIA,NNLEV)
      INTEGER MATDAT(SZN2MI,NNLEV)
      INTEGER VECDAT(SZN2VI,NNLEV)

      INTEGER IPARAM(SZISDI)
      DOUBLE PRECISION DPARAM(SZISDD)
      
      INTEGER ISTPAR(SZNSDI)
      DOUBLE PRECISION DSTPAR(SZNSDD)
      
      INTEGER IPPDAT(SZISPI)
      DOUBLE PRECISION DPPDAT(SZISPD)
      
C     Multigrid structures provide space for one Solver, 
C     MG-Solver/preconditioner, Smoother and coarse grid solver:
      
      INTEGER IMGPAR (SZ020I+3*SZSLVI)
      DOUBLE PRECISION DMGPAR (SZ020D+3*SZSLVD)
      
      INTEGER IASMBL(SZASMI)
      DOUBLE PRECISION DASMBL(SZASMD)
      
      INTEGER IGEOM(*)
      DOUBLE PRECISION DGEOM(*)
      
      INTEGER IADTS(SZADTI)
      DOUBLE PRECISION DADTS(SZADTD)
      
      INTEGER ITRACK(SZTRKI)
      
      INTEGER LUP,LRHS
      
      DOUBLE PRECISION TIMENS
      EXTERNAL GNRHSV
      
      EXTERNAL I000,DFPSTD,DFGGRI,FBPNTS
      
      INTEGER I,NUVP
      DOUBLE PRECISION WCTIM,WCTINI,WCTTOT
      
C     Store the wall clock time into WCTIM. Is used to determine the
C     total time at the end of the algorithm. 

      WCTTOT = 0D0
      CALL GTMAUX (WCTIM,WCTTOT,1,0)

C     Initialize triangulations:
C     Correct NLMIN/NLMAX into bounds 1..9 if they are out of bounds.

      CALL INMTRI (MSHOW,TRIAS,NLMIN,NLMAX,0,CMESH)

C     Initialize parameter blocks of nonstationary solver
      
      CALL INISSL (MSHOW,NLMIN,NLMAX,TRIAS,
     *             MATDAT,VECDAT,IPARAM,DPARAM,
     *             ISTPAR,DSTPAR,IMGPAR,DMGPAR,
     *             IASMBL,DASMBL,IADTS,DADTS,
     *             LUP,LRHS,NUVP)
     
C     Initialise parameter block for the postprocessing routine.
C     There is no "release" routine since there is no dynamic
C     information saved there.
C     Transfer COMMON block variables of the DAT file into the
C     parameter blocks.

      CALL ININSP (IPPDAT,DPPDAT,1)
     
C     Inform the geometry routines that we are working instationary

      IGEOM(OINSGEO) = 1

C     Open the output files for tracking of forces,...

      CALL INITRK (ITRACK,-1)
      
C     Initialize matrices

      CALL GNMATD (NLMIN,NLMAX,TRIAS,
     *             MATDAT,ISTPAR,DSTPAR,IASMBL,DASMBL)
     
C     In the instationary case, geometry/boundary information are
C     maintained by the solver internally; we mustn't implement
C     them manually into the TRIAS-structures.
C
C     Generate solution vector data. Implement boundary conditions.

      CALL GNISVC (NLMIN,NLMAX,TRIAS,
     *             VECDAT,IPARAM,DPARAM,ISTPAR,DSTPAR,IASMBL,DASMBL,
     *             IGEOM,DGEOM,
     *             1,0D0,1D0,LUP)

C     Read the starting vector if necessary

      IF (ISTART.NE.0) THEN
        CALL RDSTRT (NLMIN,NLMAX,VECDAT,TRIAS,DWORK(L(LUP)),
     *               IINT,IAPR,IAVPR,DPREP,
     *               ISTART,CSTART)
      END IF

C     Switch off mesh adaption if it's not activated in the mesh-adaption
C     DAT file; saves some space...

      IF (IGAMET.EQ.0) IPARAM(OIDMESH) = 0

C     Preparation completed; stop the time necessary for the
C     preprocessing:

      CALL GTMAUX (WCTIM,WCTINI,1,1)

C     Call the solver

      TIMENS = 0D0
      
      CALL NONST2 (NLMIN,NLMAX,
     *             TRIAS,MATDAT,VECDAT,
     *             IPARAM,DPARAM,ISTPAR,DSTPAR,
     *             IMGPAR,DMGPAR,IASMBL,DASMBL,IGEOM,DGEOM,
     *             IADTS,DADTS,IPPDAT,DPPDAT,ITRACK,
     *             NUVP,DWORK(L(LUP)),TIMENS,
     *             GNRHSV,DFGGRI,DFPSTD,FBPNTS)

C     Write final solution vector to disc?

      IF (ISOL.NE.0) THEN
        I=0
        IF (ISOL.EQ.1) I=-1
        IF (ISOL.EQ.2) I=1
        CALL PPWRVC (0,VECDAT(ONEQV,NLMAX),DWORK(L(LUP)),
     *               -1,I,CSOL)
      END IF

C     Close files with results of tracking of values (drag/lift/...) 
C     over time

      CALL DONTRK (ITRACK)
      
C     Release memory, finish

      CALL DONNSP (IPPDAT)
      
      CALL DNISSL (NLMIN,NLMAX,
     *             MATDAT,VECDAT,
     *             IPARAM,DPARAM,
     *             ISTPAR,DSTPAR,
     *             LUP,LRHS)
      
      CALL DNMTRI (NLMIN,NLMAX,TRIAS)

C     Stop the complete time including the initialisation
C     and postprocessing in WCTTOT.

      CALL GTMAUX (WCTIM,WCTTOT,1,1)

C     Print statistical data / timing results

      CALL STATIS (MSHOW,WCTINI,WCTTOT,IPARAM,DPARAM)
      
      END
     