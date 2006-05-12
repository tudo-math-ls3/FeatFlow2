
      PROGRAM CC2D_MOVBC  

      IMPLICIT NONE

      INCLUDE 'cmem.inc'
      
      INCLUDE 'ctimediscr.inc'
      
      INCLUDE 'sinigeometry.inc'
      INCLUDE 'stria.inc'
      INCLUDE 'cbasicmg.inc'
      
      INCLUDE 'cparametrization.inc'
      INCLUDE 'cdiscretization.inc'
      
      INTEGER MSHOW,MDATA
      
      INTEGER LIGEOM,LDGEOM,LTRIAS
      
      CHARACTER COPTDT*60
      
C     At the very beginning initialise the pseudodynamic memory 
C     management. Nothing must be done before that!!!

      CALL ZINIT(NNWORK,'feat.msg','data/CC2D.ERR','data/CC2D.PRT',
     *                             'data/CC2D.SYS','data/CC2D.TRC') 

C     Initialize the output channels - also at the very beginning
C     before anything is done!

      CALL OINIT (MSHOW)

C     Allocate some memory for the geometry and triangulation

      CALL ZNEW(SZGEOI,3,LIGEOM , 'IGEOM ')
      CALL ZNEW(SZGEOD,1,LDGEOM , 'DGEOM ')
      CALL ZNEW(SZTRIA*NNLEV,3,LTRIAS , 'TRIAS ')

C     Read DAT files, initialize general data, create parametrization.

      CALL INMAIN (MSHOW,KWORK(L(LIGEOM)),DWORK(L(LDGEOM)))

C     Initialize triangulations.
C     Correct NLMIN/NLMAX into bounds 1..9 if they are out of bounds.
C     Remember, NLMIN/NLMAX are COMMON-block variables!
      
      CALL INMTRI (MSHOW,KWORK(L(LTRIAS)),NLMIN,NLMAX,0,CMESH)
      
C     Read parameters of the optimizer DAT file
C     Use a string variable as some compilers have problems using a
C     direct string as parameter!

      COPTDT = 'data/optimization.dat'
      MDATA = 79
      CALL RDOPT (MDATA,MSHOW,COPTDT)
      
C     Call the optimizer for the stationary problem

      CALL CC2DOS (MSHOW,LTRIAS,LIGEOM,LDGEOM)
      
C     Clean up triangulations

      CALL DNMTRI (NLMIN,NLMAX,KWORK(L(LTRIAS)))
      
C     Clean up geometry, general cleanup

      CALL DNMAIN (KWORK(L(LIGEOM)),DWORK(L(LDGEOM)))

C     Release memory of the geometry and triangulation structure, finish
      
      CALL ZDISP(0,LTRIAS , 'TRIAS ')
      CALL ZDISP(0,LIGEOM , 'IGEOM ')
      CALL ZDISP(0,LDGEOM , 'DGEOM ')
      
      END       

************************************************************************
* CC2D - start stationary solver for standard problem
*
* This routine is called by the optimization subroutine to perform
* a standard stationary calculation without any optimization.
*
* In: 
*   MSHOW  - Output level during initialization
*   HLIST  - List of handles to the parameter blocks of the different
*            solver-components
*            = THList structure
*
* Out:
*   Output data is written to files.
************************************************************************

      SUBROUTINE CC2DST (MSHOW,HLIST)
      
      IMPLICIT NONE

      INCLUDE 'cmem.inc'
      INCLUDE 'cout.inc'

      INCLUDE 'stria.inc'
      INCLUDE 'cbasicmg.inc'
      
      INCLUDE 'dstrings.inc'
      
      INCLUDE 'stracking.inc'
      
      INCLUDE 'soptpar.inc'
      
      INCLUDE 'sgridadaptmgstaticparams.inc'
      INCLUDE 'cgridadapt.inc'
      INCLUDE 'cparametrization.inc'
      INCLUDE 'cdiscretization.inc'
      INCLUDE 'clinsol_cc2d.inc'
      INCLUDE 'cpostproc.inc'
      
      INTEGER MSHOW,HLIST(SZHLST)
      
C     local variables

      INTEGER LTRIAS ,LMATDAT,LVECDAT,LISTPAR,LDSTPAR,LIMGPAR,LDMGPAR
      INTEGER LIASMBL,LDASMBL,LIGEOM ,LDGEOM ,LUP    ,LRHS   ,LAUX  
      INTEGER IBDTYP ,IFXOUT  
      INTEGER ITRACK(SZTRKI),LFNAMS(6)
      INTEGER NUVP
      
      INTEGER I,KTRIA
      DOUBLE PRECISION WCTIM,WCTTOT,WCTINI

C     Store the wall clock time into WCTIM. Is used to determine the
C     total time at the end of the algorithm. 

      WCTTOT = 0D0
      CALL GTMAUX (WCTIM,WCTTOT,1,0)

C     Initialize everything which is left
C     for the solution process.
C
C     At first fetch the handles to local variables

      LTRIAS  = HLIST(OLTRIAS )
      LMATDAT = HLIST(OLMATDAT)
      LVECDAT = HLIST(OLVECDAT)
      LISTPAR = HLIST(OLISTPAR)
      LDSTPAR = HLIST(OLDSTPAR)
      LIMGPAR = HLIST(OLIMGPAR)
      LDMGPAR = HLIST(OLDMGPAR)
      LIASMBL = HLIST(OLIASMBL)
      LDASMBL = HLIST(OLDASMBL)
      LIGEOM  = HLIST(OLIGEOM )
      LDGEOM  = HLIST(OLDGEOM )
      LUP     = HLIST(ODUP    )
      LRHS    = HLIST(ODRHS   )
      
      NUVP    = HLIST(ONUVP  ) 
      
      IFXOUT  = HLIST(OIFXOUT)

C     Adapt the grid if activated in the user

      IF (IGAMET.EQ.3) THEN
        CALL XSMMGW (KWORK(L(LTRIAS)),NLMIN,NLMAX,
     *               KWORK(L(LIGEOM)),DWORK(L(LDGEOM)),
     *               0,0D0)
      END IF
      
C     Initialize/Generate matrices
     
      CALL GNMATD (NLMIN,NLMAX,KWORK(L(LTRIAS)),
     *             KWORK(L(LMATDAT)),KWORK(L(LISTPAR)),
     *             DWORK(L(LDSTPAR)),KWORK(L(LIASMBL)),
     *             DWORK(L(LDASMBL)))
     
C     Generate boundary/geometry information for the discretization.
C     Generate precalculated infomration about geometries.
C     As we are in the stationary case, we implement boundary 
C     conditions directly into the triangulation.
     
      CALL GNBDGE (NLMIN,NLMAX,KWORK(L(LTRIAS)),
     *             KWORK(L(LIASMBL)),DWORK(L(LDASMBL)),
     *             KWORK(L(LIGEOM)),DWORK(L(LDGEOM)))

      IBDTYP = 0
      IF (IBDR.EQ.1) IBDTYP = 1
      IBDTYP = IOR(4+2,IBDTYP)
     
C     Generate RHS/solution vector data
     
      CALL GNSTVC (NLMIN,NLMAX,KWORK(L(LTRIAS)),
     *             KWORK(L(LVECDAT)),KWORK(L(LISTPAR)),
     *             DWORK(L(LDSTPAR)),KWORK(L(LIASMBL)),
     *             DWORK(L(LDASMBL)),
     *             KWORK(L(LIGEOM)),DWORK(L(LDGEOM)),
     *             IBDTYP,LUP,LRHS)

C     Read the starting vector if necessary

      IF (ISTART.NE.0) THEN
        CALL RDSTRT (NLMIN,NLMAX,KWORK(L(LVECDAT)),KWORK(L(LTRIAS)),
     *               DWORK(L(LUP)),
     *               IINT,IAPR,IAVPR,DPREP,
     *               ISTART,CSTART)
      END IF

C     Preparation completed; stop the time necessary for the
C     preprocessing:

      CALL GTMAUX (WCTIM,WCTINI,1,1)

C     Call the solver

      CALL NSDEF2 (NLMIN,NLMAX,KWORK(L(LTRIAS)),
     *             KWORK(L(LMATDAT)),KWORK(L(LVECDAT)),
     *             KWORK(L(LISTPAR)),DWORK(L(LDSTPAR)),
     *             KWORK(L(LIMGPAR)),DWORK(L(LDMGPAR)),
     *             KWORK(L(LIASMBL)),DWORK(L(LDASMBL)),
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
     *             KWORK(L(LTRIAS)),KWORK(L(LMATDAT)),KWORK(L(LVECDAT)),
     *             KWORK(L(LISTPAR)),DWORK(L(LDSTPAR)),
     *             KWORK(L(LIMGPAR)),DWORK(L(LDMGPAR)),
     *             KWORK(L(LIASMBL)),DWORK(L(LDASMBL)),
     *             KWORK(L(LIGEOM)),DWORK(L(LDGEOM)),
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
     *             KWORK(L(LTRIAS)),KWORK(L(LIASMBL)),DWORK(L(LDASMBL)))
      
C     Stop the complete time including the initialisation
C     and postprocessing in WCTTOT.

      CALL GTMAUX (WCTIM,WCTTOT,1,1)

C     Print statistical data / timing results

      CALL STATST (MSHOW,WCTINI,WCTTOT,
     *             KWORK(L(LISTPAR)),DWORK(L(LDSTPAR)))
      
1     FORMAT(79('-'))

      END
      
************************************************************************
* Probe stationary solver
*
* This routine is a callback routine for all optimization algorithms.
* It calls the stationary solver to compute a quantity. The desired
* quantity is returned in the parameters.
*
* The routine assumes IINFO to be a THList structure filled with
* data about the problem, parameters, triangulations and so on.
*
* In: 
*   IDIMS   - Number of dimensions of the design parameters
*   LCOORD  - Handle to array [1..IDIMS] of double
*             Design parameter setting where to evaluate
*   NRES    - integer
*             Size of result array
*   LRES    - integer
*             Handle to array [1..NRES] of double,
*             where to write the results to
*   PITER   - Number of current probing iterate
*   IINFO   - THList structure with information how to solve
*             the problem
*   DINFO   - Double precision parameter block.
*             Not used here.
*     
*
* Out:
*   RES     - array [1..NRES] of double
*             RES(1) is the value of the point DCOORDS,
*             which is used by the optimizer to
*             optimize for. RES(2..NRES) can be set by
*             DPROBE to arbitrary values that are later
*             passed to DNXSTP as result values.
************************************************************************

      SUBROUTINE DPRSTS (IDIMS,LCOORDS,NRES,LRES,PITER,IINFO,DINFO)
      
      IMPLICIT NONE

      INCLUDE 'cmem.inc'
      INCLUDE 'cout.inc'

      INCLUDE 'stria.inc'
      INCLUDE 'cbasicmg.inc'
      
      INCLUDE 'dstrings.inc'
      
      INCLUDE 'sassembly.inc'
      INCLUDE 'sinigeometry.inc'
      
      INCLUDE 'soptpar.inc'
      INCLUDE 'sinioptimization.inc'
      
      INCLUDE 'sgridadaptmgstaticparams.inc'
      INCLUDE 'cgridadapt.inc'
      INCLUDE 'cparametrization.inc'
      INCLUDE 'cdiscretization.inc'
      INCLUDE 'clinsol_cc2d.inc'
      INCLUDE 'cpostproc.inc'
      
      INCLUDE 'sgeometries.inc'
      INCLUDE 'scompgeometries.inc'
      
      INTEGER IDIMS,LCOORDS,NRES,LRES,PITER,IINFO(*)
      DOUBLE PRECISION DINFO(*)
      
C     local variables

      INTEGER LTRIAS ,LMATDAT,LVECDAT,LISTPAR,LDSTPAR,LIMGPAR,LDMGPAR
      INTEGER LIASMBL,LDASMBL,LIGEOM ,LDGEOM ,LUP    ,LRHS   
      INTEGER IBDTYP ,IFXOUT,LIOPTP,LDOPTP
      INTEGER LCURTR ,LAUX ,I,IOPTTP ,IOPTPR, L1, L2, KIGEOM,KDGEOM
      DOUBLE PRECISION COEF1,COEF2,RE1
      DOUBLE PRECISION R1,R2,RSTMP(10)
      INTEGER NUVP
      
      DOUBLE PRECISION PI
      PARAMETER (PI=3.1415926535897931D0)
      
C     Initialize everything which is left
C     for the solution process.
C
C     At first fetch the handles to local variables

      LTRIAS  = IINFO(OLTRIAS )
      LCURTR  = IINFO(OLCURTR )
      LMATDAT = IINFO(OLMATDAT)
      LVECDAT = IINFO(OLVECDAT)
      LISTPAR = IINFO(OLISTPAR)
      LDSTPAR = IINFO(OLDSTPAR)
      LIMGPAR = IINFO(OLIMGPAR)
      LDMGPAR = IINFO(OLDMGPAR)
      LIASMBL = IINFO(OLIASMBL)
      LDASMBL = IINFO(OLDASMBL)
      LIGEOM  = IINFO(OLIGEOM )
      LDGEOM  = IINFO(OLDGEOM )
      LUP     = IINFO(ODUP    )
      LRHS    = IINFO(ODRHS   )
      
      NUVP    = IINFO(ONUVP  )
      
      IFXOUT  = IINFO(OIFXOUT)
      
      LIOPTP  = IINFO(OLIOPTP)
      LDOPTP  = IINFO(OLDOPTP)
      IOPTTP  = KWORK(L(LIOPTP)+OIOPTTP-1)
      IOPTPR  = KWORK(L(LIOPTP)+OIOPTPR-1)
      
      KIGEOM  = L(LIGEOM)
      KDGEOM  = L(LDGEOM)
      
      IF (IOPTPR.EQ.0) THEN
      
C       Adapt the geometry to the current situation.
C       In LCOORDS we can find the new X-, Y-position and radius
C       of the first fictitious boundary object. Store this
C       info in the geometry structure.
C       The array of LCOORDS is set up in CC2DOS!
C
C       The position of the object is depending on the type of
C       simulation. For IOPTTP=2, we have a standard brute-force
C       test with all parameter settings.
C       For IOPTTP >= 100, we are using a real optimization algorithm
C       which adapts the position of the object.

        IF (IOPTTP.EQ.2) THEN

          DWORK(L(LDGEOM)+ODCXPOS-1) = DWORK(L(LCOORDS))
          DWORK(L(LDGEOM)+ODCYPOS-1) = DWORK(L(LCOORDS)+1)
          DWORK(L(LDGEOM)+ODCROT -1) = DWORK(L(LCOORDS)+2)
          DWORK(L(LDGEOM)+ODCRSIN-1) = SIN(DWORK(L(LCOORDS)+2)*PI/180.0)
          DWORK(L(LDGEOM)+ODCRCOS-1) = COS(DWORK(L(LCOORDS)+2)*PI/180.0)
        
        ELSE IF (IOPTTP.GE.100) THEN
        
C         For IOPTTP >= 100, we are using a real optimization algorithm
C         which adapts the position of the object. This is the same
C         formula as in the case IOPTTP=2:

          DWORK(L(LDGEOM)+ODCXPOS-1) = DWORK(L(LCOORDS))
          DWORK(L(LDGEOM)+ODCYPOS-1) = DWORK(L(LCOORDS)+1)
          DWORK(L(LDGEOM)+ODCROT -1) = DWORK(L(LCOORDS)+2)
          DWORK(L(LDGEOM)+ODCRSIN-1) = SIN(DWORK(L(LCOORDS)+2)*PI/180.0)
          DWORK(L(LDGEOM)+ODCRCOS-1) = COS(DWORK(L(LCOORDS)+2)*PI/180.0)
        
        END IF
        
        CALL UPDGEO (DWORK(L(LIGEOM)),DWORK(L(LDGEOM)))
        
C       If the geometry parameters are out of range, cancel!

        IF ( (DWORK(L(LCOORDS))  .GT.DWORK(L(LDOPTP)+OOBXMAX-1)) .OR.
     *       (DWORK(L(LCOORDS)+1).GT.DWORK(L(LDOPTP)+OOBYMAX-1)) .OR.
     *       (DWORK(L(LCOORDS)+2).GT.DWORK(L(LDOPTP)+OOBRMAX-1)) .OR.
     *       (DWORK(L(LCOORDS))  .LT.DWORK(L(LDOPTP)+OOBXMIN-1)) .OR.
     *       (DWORK(L(LCOORDS)+1).LT.DWORK(L(LDOPTP)+OOBYMIN-1)) .OR.
     *       (DWORK(L(LCOORDS)+2).LT.DWORK(L(LDOPTP)+OOBRMIN-1)) ) THEN
          DWORK(L(LRES)) = 1D99
          RETURN
        END IF
      
      ELSE IF (IOPTPR.EQ.10) THEN
      
C       This is a slightly more complex situation:
C       We have two balls and must correctly set their Y-positions!
C       DCOORDS(1) and DCOORDS(2) define the both Y-positions.
C
C       Remember the starting address of OIGEODL in DWORK.

        L1 = KWORK(KIGEOM+OIGEODL-1)
        IF (L1.EQ.0) THEN
          WRITE (*,*) 'Geometry data array IGEODL undefined!'
          STOP
        END IF
        L1 = L(L1)

C       Seek to the structure of the first circle

        L2 = KWORK(KIGEOM+OIGEOMS-1+3*0+1)
        L2 = L1+L2-1
        
C       Update the coordinate system of the object

        CALL GUPCOR (DWORK(L2),
     *       DWORK(KDGEOM+ODCXPOS-1), DWORK(L(LCOORDS)), 0D0, 1D0)

C       What's the radius of the circle?

        R1 = DWORK(L2+OCRAD-1)

C       Seek to the structure of the 2nd circle

        L2 = KWORK(KIGEOM+OIGEOMS-1+3*1+1)
        L2 = L1+L2-1
        
C       Update its coordinate system

        CALL GUPCOR (DWORK(L2),
     *       DWORK(KDGEOM+ODCXPOS-1), DWORK(L(LCOORDS)+1), 0D0, 1D0)

C       What's the radius of the circle?

        R2 = DWORK(L2+OCRAD-1)
      
C       If the geometry parameters are out of range
C       or if the circles are overlapping, cancel!
      
        IF ( ABS(DWORK(L(LCOORDS))-DWORK(L(LCOORDS)+1))
     *       .LT.(R1+R2) ) THEN
          DWORK(L(LRES)) = 1D99 
          RETURN
        END IF
      
        IF ( (DWORK(L(LCOORDS))  .GT.DWORK(L(LDOPTP)+OOBYMAX-1)-R1) .OR.
     *       (DWORK(L(LCOORDS)+1).GT.DWORK(L(LDOPTP)+OOBYMAX-1)-R2) .OR.
     *       (DWORK(L(LCOORDS))  .LT.DWORK(L(LDOPTP)+OOBYMIN-1)+R1) .OR.
     *       (DWORK(L(LCOORDS)+1).LT.DWORK(L(LDOPTP)+OOBYMIN-1)+R2) ) 
     *  THEN
          DWORK(L(LRES)) = 1D99
          RETURN
        END IF
      
      END IF

C     Adapt the grid if activated in the user

      IF (IGAMET.EQ.3) THEN
      
C       Restore original grid.
C       Go through all levels and call the grid-restore routine
C       to restore the original mesh.
C       Restore everything we duplicated.

        DO I=NLMIN,NLMAX
          CALL TRIRST (KWORK(L(LCURTR)+SZTRIA*(I-1)),
     *                 KWORK(L(LTRIAS)+SZTRIA*(I-1)))
        END DO
      
C       Adapt the grid in "CURTRI"
      
        CALL XSMMGW (KWORK(L(LCURTR)),NLMIN,NLMAX,
     *               KWORK(L(LIGEOM)),DWORK(L(LDGEOM)))
     
      END IF
      
C     Initialize/Generate matrices
     
      CALL GNMATD (NLMIN,NLMAX,KWORK(L(LCURTR)),
     *             KWORK(L(LMATDAT)),KWORK(L(LISTPAR)),
     *             DWORK(L(LDSTPAR)),KWORK(L(LIASMBL)),
     *             DWORK(L(LDASMBL)))
     
C     Generate boundary/geometry information for the discretization.
C     As we are in the stationary case, we implement boundary 
C     conditions directly into the triangulation.
     
      CALL GNBDGE (NLMIN,NLMAX,KWORK(L(LCURTR)),
     *             KWORK(L(LIASMBL)),DWORK(L(LDASMBL)),
     *             KWORK(L(LIGEOM)),DWORK(L(LDGEOM)))

      IBDTYP = 0
      IF (IBDR.EQ.1) IBDTYP = 1
      IBDTYP = IOR(4+2,IBDTYP)
     
C     Read the starting vector if necessary or fill it with 0

      IF (ISTART.NE.0) THEN
        CALL RDSTRT (NLMIN,NLMAX,KWORK(L(LVECDAT)),KWORK(L(LCURTR)),
     *               DWORK(L(LUP)),
     *               IINT,IAPR,IAVPR,DPREP,
     *               ISTART,CSTART)
      ELSE
        CALL LCL1(DWORK(L(LUP)),NUVP)
      END IF

C     Generate RHS vector data.
C     Implement boundary conditions into RHS and solution.
     
      CALL GNSTVC (NLMIN,NLMAX,KWORK(L(LCURTR)),
     *             KWORK(L(LVECDAT)),KWORK(L(LISTPAR)),
     *             DWORK(L(LDSTPAR)),KWORK(L(LIASMBL)),
     *             DWORK(L(LDASMBL)),
     *             KWORK(L(LIGEOM)),DWORK(L(LDGEOM)),
     *             IBDTYP,LUP,LRHS)

C     Call the solver

      CALL NSDEF2 (NLMIN,NLMAX,KWORK(L(LCURTR)),
     *             KWORK(L(LMATDAT)),KWORK(L(LVECDAT)),
     *             KWORK(L(LISTPAR)),DWORK(L(LDSTPAR)),
     *             KWORK(L(LIMGPAR)),DWORK(L(LDMGPAR)),
     *             KWORK(L(LIASMBL)),DWORK(L(LDASMBL)),
     *             NUVP,
     *             DWORK(L(LUP)),DWORK(L(LRHS)))

      IF (IOPTPR.EQ.0) THEN

C       We want to calculate forces on the first occurring fictitious
C       boundary object. 
C       Get the coefficients of the integral:

        RE1 = DWORK(L(LDASMBL)+ORE-1)
        CALL FPTSIN(5,0,1D0,RE1,
     *              KWORK(L(LIASMBL)),DWORK(L(LDASMBL)),
     *              KWORK(L(LIGEOM)),DWORK(L(LDGEOM)),
     *              KWORK(L(LCURTR)+SZTRIA*(NLMAX-1)),
     *              0,COEF1,COEF2)
        
C       Call the calculation method. Store the values to RES(2..*)

        CALL DFPVIN (NLMIN,NLMAX,
     *               KWORK(L(LCURTR)),
     *               KWORK(L(LMATDAT)),KWORK(L(LVECDAT)),
     *               KWORK(L(LISTPAR)),DWORK(L(LDSTPAR)),
     *               KWORK(L(LIMGPAR)),DWORK(L(LDMGPAR)),
     *               KWORK(L(LIASMBL)),DWORK(L(LDASMBL)),
     *               KWORK(L(LIGEOM)),DWORK(L(LDGEOM)),
     *               NUVP,DWORK(L(LUP)),DWORK(L(LRHS)),
     *               1, 0, 0D0,0D0, 4, COEF1,COEF2,
     *               DWORK(L(LRES)+1))
     
C       Put the value of interest to RES(1) - we take the absolute lift

        DWORK(L(LRES)) = ABS(DWORK(L(LRES)+2))
        
      ELSE IF (IOPTPR.EQ.10) THEN
      
C       Calculate the absolute lift-forces of both objects.
C       Don't calculate the coefficients, directly calculate the
C       forces (-> set dpf1=1, dpf2=2 to directly calculate forces).

        CALL DFPVIN (NLMIN,NLMAX,
     *               KWORK(L(LCURTR)),
     *               KWORK(L(LMATDAT)),KWORK(L(LVECDAT)),
     *               KWORK(L(LISTPAR)),DWORK(L(LDSTPAR)),
     *               KWORK(L(LIMGPAR)),DWORK(L(LDMGPAR)),
     *               KWORK(L(LIASMBL)),DWORK(L(LDASMBL)),
     *               KWORK(L(LIGEOM)),DWORK(L(LDGEOM)),
     *               NUVP,DWORK(L(LUP)),DWORK(L(LRHS)),
     *               1, 2, 0D0,0D0, 4, 1D0,2D0,
     *               RSTMP)
     
C       Add a little bit gravitation to reduce the lift...
     
        DWORK(L(LRES)+1) = RSTMP(2)
C         - 20.0 * 9.81*PI*R1**2

        CALL DFPVIN (NLMIN,NLMAX,
     *               KWORK(L(LCURTR)),
     *               KWORK(L(LMATDAT)),KWORK(L(LVECDAT)),
     *               KWORK(L(LISTPAR)),DWORK(L(LDSTPAR)),
     *               KWORK(L(LIMGPAR)),DWORK(L(LDMGPAR)),
     *               KWORK(L(LIASMBL)),DWORK(L(LDASMBL)),
     *               KWORK(L(LIGEOM)),DWORK(L(LDGEOM)),
     *               NUVP,DWORK(L(LUP)),DWORK(L(LRHS)),
     *               1, 3, 0D0,0D0, 4, 1D0,2D0,
     *               RSTMP)
     
        DWORK(L(LRES)+2) = RSTMP(2)
C         - 20.0 * 9.81*PI*R2**2
        
C       The sum of both is our fitness value

C        DWORK(L(LRES)) = MAX(0.01,ABS(DWORK(L(LRES)+1))) + 
C     *                   MAX(0.01,ABS(DWORK(L(LRES)+2)))
        DWORK(L(LRES)) = ABS(DWORK(L(LRES)+1))**2 + 
     *                   ABS(DWORK(L(LRES)+2))**2
      
      END IF

C     Memory is released in the DNXSTS routine after postprocessing

1     FORMAT(79('-'))

      END
      
************************************************************************
* Cleanup stationary solver, go to next optimization step
*
* This routine is called by the optimizer after DPROBE. It informs
* whether or not the new iterate is accepted.
* The routine performs postprocessing with the results of DPROBE
* and cleanup of the structures that were allocated in DPROBE.
*
* In:
*   DCOORDS - array [1..IDIMS] of double precision
*             Parameter set of probed iterate.
*   DCOOPT  - array [1..IDIMS] of double precision
*             Parameter set of current optimal iterate
*   IINFO,
*   DINFO   - Integer / double precision parameter
*             block from above
*   PITER   - Number of current probing iterate
*   STEP    - Number of current accepted iterate
*   IACPT   - whether DCOORDS will be accepted as the new
*             iterate.
*             =0: DCOORDS will not be accepted
*             =1: DCOORDS will be accepted
*             Can be modified by DNXSTP to prevent the
*             point from being accepted.
*   NRES    - integer
*             Size of result array RES
*   DRES    - Result array that was calculated 
*             in DPROBE for the parameter setting DCOORDS.
*   DRSOPT  - Result array that was calculated 
*             in DPROBE for the parameter setting DCOOPT.
*
* Out:
*   -
************************************************************************
      
      SUBROUTINE DNXSTS (IDIMS,DCOORDS,DCOOPT,IINFO,DINFO,
     *                   PITER,STEP,IACPT,NRES,DRES,DRSOPT)
     
      IMPLICIT NONE
      
      INTEGER IDIMS,IINFO(*),PITER,STEP,IACPT,NRES
      DOUBLE PRECISION DCOORDS(IDIMS),DCOOPT(IDIMS)
      DOUBLE PRECISION DINFO(*),DRES(NRES),DRSOPT(NRES)
      
      INCLUDE 'cmem.inc'
      INCLUDE 'cerr.inc'
      INCLUDE 'cout.inc'

      INCLUDE 'stria.inc'
      INCLUDE 'cbasicmg.inc'
      
      INCLUDE 'dstrings.inc'
      
      INCLUDE 'stracking.inc'
      INCLUDE 'sassembly.inc'
      INCLUDE 'sinioptimization.inc'
      INCLUDE 'sinigeometry.inc'
      
      INCLUDE 'soptpar.inc'
      
      INCLUDE 'sgridadaptmgstaticparams.inc'
      INCLUDE 'cgridadapt.inc'
      INCLUDE 'cparametrization.inc'
      INCLUDE 'cdiscretization.inc'
      INCLUDE 'clinsol_cc2d.inc'
      INCLUDE 'cpostproc.inc'     
      
      INCLUDE 'sgeometries.inc'
      
      INCLUDE 'slinsol.inc'
      INCLUDE 'slinsol2dns.inc'
      
C     local variables

      INTEGER LTRIAS ,LMATDAT,LVECDAT,LISTPAR,LDSTPAR,LIMGPAR,LDMGPAR
      INTEGER LIASMBL,LDASMBL,LIGEOM ,LDGEOM ,LUP    ,LRHS   ,LAUX  
      INTEGER IFXOUT ,LCURTR ,LIOPTP ,IOPTTP ,IOPTPR
      INTEGER ITRACK(SZTRKI),LFNAMS(6),KIGEOM
      INTEGER NUVP,NEQU
      
      INTEGER I,K,L2,L1
      DOUBLE PRECISION Y1,Y2
      
      EXTERNAL DFGMMV,DFGMMC,UE
      DOUBLE PRECISION UE

      CHARACTER CFN*60

C     Initialize everything which is left
C     for the solution process.
C
C     At first fetch the handles to local variables

      LTRIAS  = IINFO(OLTRIAS )
      LCURTR  = IINFO(OLCURTR )
      LMATDAT = IINFO(OLMATDAT)
      LVECDAT = IINFO(OLVECDAT)
      LISTPAR = IINFO(OLISTPAR)
      LDSTPAR = IINFO(OLDSTPAR)
      LIMGPAR = IINFO(OLIMGPAR)
      LDMGPAR = IINFO(OLDMGPAR)
      LIASMBL = IINFO(OLIASMBL)
      LDASMBL = IINFO(OLDASMBL)
      LIGEOM  = IINFO(OLIGEOM )
      LDGEOM  = IINFO(OLDGEOM )
      LUP     = IINFO(ODUP    )
      LRHS    = IINFO(ODRHS   )
      
      NUVP    = IINFO(ONUVP  )
      
      IFXOUT  = IINFO(OIFXOUT)
      
      LIOPTP  = IINFO(OLIOPTP)
      IOPTTP  = KWORK(L(LIOPTP)+OIOPTTP-1)
      IOPTPR  = KWORK(L(LIOPTP)+OIOPTPR-1)
      
      KIGEOM  = L(LIGEOM)
      
      NEQU    = KWORK(L(LVECDAT)+(NLMAX-1)*SZN2VI+ONU-1)
      
C     In the very first iteration, write the headline of 
C     the output file

      IF (IOPTPR.EQ.0) THEN
        
C       Standard problem: one ball in a channel, determine drag/lift)
        
        IF (PITER.EQ.1) THEN
        
          WRITE (IINFO(OIFXOUT),'(A)') 
     *                   '#     X-Pos                   '//
     *                   '  Y-Pos                   '//
     *                   '  Rotation               '//
     *                   'Drag                    '//
     *                   'Lift                    '//
     *                   'abs.Lift                '//
     *                   'Accepted'
     
        END IF
     
C       Print the results of the previous computation:

        WRITE (MTERM,'(A)')      '-------------------------------------'
        WRITE (MTERM,'(A,I10)')    'Iteration           : ',PITER
        WRITE (MTERM,'(A,I10)')    'Acceptance counter  : ',STEP
        WRITE (MTERM,'(A,3F13.6)') 'Position            : ',
     *                             DCOORDS(1),DCOORDS(2),DCOORDS(3)
        WRITE (MTERM,'(A,D24.10)') 'Probing value       : ',DRES(1)
        WRITE (MTERM,'(A,D24.10)') 'Best probing value  : ',DRSOPT(1)
        WRITE (MTERM,'(A,I10)')    'Accepted            : ',IACPT
        WRITE (MTERM,'(A)')      '-------------------------------------'

C       Cancel if the value is 1D99 - there's nothing calculated!
        
        IF (DRES(1).EQ.1D99) RETURN

C       Write out the results for drag/lift to a file

        WRITE (IFXOUT,'(F24.16,2X,F24.16,2X,F24.16,3E24.13,I6)') 
     *            DWORK(L(LDGEOM)+ODCXPOS-1),
     *            DWORK(L(LDGEOM)+ODCYPOS-1),
     *            DWORK(L(LDGEOM)+ODCROT-1),
     *            REAL(DRES(2)),REAL(DRES(3)),REAL(ABS(DRES(3))),IACPT
        CALL FLUSH(IFXOUT)

      ELSE IF (IOPTPR.EQ.10) THEN
        
C       Two balls in a channel, determine sum of abs. lifts

        IF (PITER.EQ.1) THEN
        
          WRITE (IINFO(OIFXOUT),'(A)') 
     *                   '#     Y-Pos_1                 '//
     *                   '  Y-Pos_2                 '//
     *                   'Lift_1                   '//
     *                   'Lift_2                  '//
     *                   'Sum(abs.Lifts)          '//
     *                   'Accepted'
     
        END IF
      
C       Print the results of the previous computation:

        WRITE (MTERM,'(A)')      '-------------------------------------'
        WRITE (MTERM,'(A,I10)')    'Iteration           : ',PITER
        WRITE (MTERM,'(A,I10)')    'Acceptance counter  : ',STEP
        WRITE (MTERM,'(A,2F13.6)') 'Position            : ',
     *                             DCOORDS(1),DCOORDS(2)
        WRITE (MTERM,'(A,3F13.6)') 'Lift:               : ',
     *                             DRES(2),DRES(3)
        WRITE (MTERM,'(A,D24.10)') 'Probing value       : ',DRES(1)
        WRITE (MTERM,'(A,D24.10)') 'Best probing value  : ',DRSOPT(1)
        WRITE (MTERM,'(A,I10)')    'Accepted            : ',IACPT
        WRITE (MTERM,'(A)')      '-------------------------------------'

C       Cancel if the value is 1D99 - there's nothing calculated!
        
        IF (DRES(1).EQ.1D99) RETURN

        L1 = L(KWORK(KIGEOM+OIGEODL-1))

        L2 = KWORK(KIGEOM+OIGEOMS-1+3*0+1)
        Y1 = DWORK(L1+L2-1+OREFY-1)
        L2 = KWORK(KIGEOM+OIGEOMS-1+3*1+1)
        Y2 = DWORK(L1+L2-1+OREFY-1)

        WRITE (IFXOUT,'(F24.16,2X,F24.16,2X,3E24.13,I6)') 
     *            Y1,Y2,
     *            REAL(DRES(2)),REAL(DRES(3)),REAL(DRES(1)),IACPT
        CALL FLUSH(IFXOUT)

      END IF

C     Call the postprocessing if the step is accepted or if we
C     perform brute-force tests.
C     IOPTTP=1..49 are reserved for brute force tests of different
C     type.

      IF ((IACPT.NE.0).OR.(IINFO(OIOPTTP).LT.100)) THEN

C       Postprocessing:
C
C       Open tracking files
      
        CALL INITRK (ITRACK,-1)
      
        WRITE (MTERM,'(A)') 'Starting postprocessing...'

C       Reserve memory for an auxiliary vector beforehand.

        CALL ZNEW (NUVP,-1,LAUX,'DAUX  ')
        
C       Prepare file output
        
        I=0
        IF (ISOL.EQ.1) I=-1
        IF (ISOL.EQ.2) I=1
        
C       Fetch names of output files into DStrings

        LFNAMS(1) = STNEWC (.TRUE.,CSOL)
        LFNAMS(2) = STNEWC (.TRUE.,CFLDU )
        LFNAMS(3) = STNEWC (.TRUE.,CFLDV )
        LFNAMS(4) = STNEWC (.TRUE.,CFLDP )
        LFNAMS(5) = STNEWC (.TRUE.,CFLISO)
        LFNAMS(6) = STNEWC (.TRUE.,CFLGMV)
        
C       Call the postprocessing routine. Don't write GMV-files!
        
        CALL DFPSTA (NLMIN,NLMAX,
     *               KWORK(L(LCURTR)),KWORK(L(LMATDAT)),
     *               KWORK(L(LVECDAT)),
     *               KWORK(L(LISTPAR)),DWORK(L(LDSTPAR)),
     *               KWORK(L(LIMGPAR)),DWORK(L(LDMGPAR)),
     *               KWORK(L(LIASMBL)),DWORK(L(LDASMBL)),
     *               KWORK(L(LIGEOM)),DWORK(L(LDGEOM)),
     *               ITRACK,NUVP,DWORK(L(LUP)),
     *               DWORK(L(LRHS)),DWORK(L(LAUX)),
     *               0,IERANA,1,4,1,
     *               I*NLMAX,IFUSAV,IFPSAV,IFXSAV,LFNAMS)

C       We write GMV-files by hand.

        IF (IGMV.GT.0) THEN
        
C         We want to write the GMV file to disc on level IGMV, i.e. at
C         level:

          I = MAX(NLMIN,MIN(NLMAX,IGMV))

C         -------------------------------------------------------------
C         Use handle 69 for writing the GMV file.

          K = 69
          
C         Create a name for the GMV file.

          CALL STPUT (LFNAMS(6), CFN)
          
C         In brute force tests we name the file for the iteration
C         counter, otherwise for the number of the accepted step.
          
          IF (IOPTTP.GE.100) THEN
            CALL GMVOF0 (K,STEP,CFN)
          ELSE
            CALL GMVOF0 (K,PITER,CFN)
          END IF
          
          IF (IER.EQ.0) THEN
          
C           Write a standard GMV file: Header, content, footer.

            CALL GMVHEA (K)
            CALL AUXGMV (K,KWORK(L(LCURTR)),NLMAX,I,NEQU,DWORK(L(LUP)),
     *                   DWORK(L(LDASMBL)+ORE-1),
     *                   KWORK(L(LIASMBL)),DWORK(L(LDASMBL)),
     *                   KWORK(L(LIGEOM)),DWORK(L(LDGEOM)))
            CALL GMVFOT (K)
          
            CLOSE (K)
          
C           -------------------------------------------------------------
C           That's it, GMV file is written out.
      
          END IF ! IER = 0
          
        END IF

C       Release memory for postprocessing
        
        DO I=6,1,-1
          CALL STDIS(LFNAMS (I))
        END DO

        CALL ZDISP (0,LAUX,'DAUX  ')
        
C       Close tracking files

        CALL DONTRK (ITRACK)     
        
      END IF

C     Release memory, finish

      CALL DNBDGE (NLMIN,NLMAX,
     *             KWORK(L(LCURTR)),KWORK(L(LIASMBL)),DWORK(L(LDASMBL)))
      
C      IACPT = 1
      
1     FORMAT(79('-'))

      END
     
      
************************************************************************
* CC2DOS - start stationary solver for standard problem
*          with optimization extension
*
* LTRIAS - Handle to array[1..SZTRIA,1..NNLEV]
*          Triangulation structures on all levels
* LIGEOM,
* LDGEOM - Handles to integer/double parameter blocks of the 
*          (moving) geometry
************************************************************************

      SUBROUTINE CC2DOS (MSHOW,LTRIAS,LIGEOM,LDGEOM)
      
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
      
      INCLUDE 'cfileout.inc'
      
      INCLUDE 'cparametrization.inc'
      INCLUDE 'cdiscretization.inc'
      INCLUDE 'sgridadaptmgstaticparams.inc'
      INCLUDE 'cgridadapt.inc'
      
      INCLUDE 'ciniopt.inc'
      INCLUDE 'sinioptimization.inc'
      
      INCLUDE 'soptpar.inc'

      INTEGER MSHOW,LIGEOM,LDGEOM,LTRIAS
      
C     We need a *couple* of structures that have to be initialized
C     and passed to the later optimizer, so the optimizer can pass
C     them to the solver for calculat the desired values.
C     The problem is how to pass them, since the optimizer "only"
C     accepts an integer and a double precision data array.
C     We can of course set up one large parameter block with all
C     sub-parameter blocks in there, but that would end up in a
C     mess.
C     We make another approach: For every parameter block we allocate
C     memory on the heap. Then, we pass all handles to the optimizer.
C     The optimizer passes them to the solver, so the solver can
C     access all parameters.
C
C     The structure of the parameter list is documented in the
C     "soptpar.inc" file.
C
      INTEGER HLIST (SZHLST)
      INTEGER NDIM
      
C     DLIST is a dummy.
      
      DOUBLE PRECISION DLIST
      
      CHARACTER CSTR*255,CBFN*60
      INTEGER IFMTS,I,J,LIOPTP,LDOPTP
      
C     Variables where the optimizer stores intermediate results

      DOUBLE PRECISION DSTART(10),DSTEP(10),DEND(10),DFINAL(10)
      DOUBLE PRECISION FINTOL(10)
      INTEGER STEPS
      
      EXTERNAL DNXSTS,DPRSTS
      
C     For easier access:

      INTEGER LMATDAT,LVECDAT,LCURTR
      INTEGER LISTPAR,LDSTPAR,LIMGPAR,LDMGPAR,LIASMBL,LDASMBL
      INTEGER LUP,LRHS,NUVP
      
C     Allocate memory for all the quantities:

      CALL ZNEW(SZTRIA*NNLEV,3,LCURTR , 'CURTRI')
      CALL ZNEW(SZN2MI*NNLEV,3,LMATDAT, 'MATDAT')
      CALL ZNEW(SZN2VI*NNLEV,3,LVECDAT, 'VECDAT')
      CALL ZNEW(SZNSDI,3,LISTPAR, 'ISTPAR')
      CALL ZNEW(SZNSDD,1,LDSTPAR, 'DSTPAR')
      CALL ZNEW(SZ020I+3*SZSLVI,3,LIMGPAR, 'IMGPAR')
      CALL ZNEW(SZ020D+3*SZSLVD,1,LDMGPAR, 'DMGPAR')
      CALL ZNEW(SZASMI,3,LIASMBL, 'IASMBL')
      CALL ZNEW(SZASMD,1,LDASMBL, 'DASMBL')
      CALL ZNEW(SZOPTI,3,LIOPTP, 'IOPTPR')
      CALL ZNEW(SZOPTD,1,LDOPTP, 'DOPTPR')

      WRITE (CSTR,9000)
      CALL CNOUTS (MSHOW,MTERM,MFILE,.TRUE.,CSTR)

C     Intialize parameter block of optimizer

      CALL INIOPT(KWORK(L(LIOPTP)),DWORK(L(LDOPTP)),1)

C     Initialize stationary solver

      CALL INSTSL (NLMIN,NLMAX,KWORK(L(LTRIAS)),
     *             KWORK(L(LMATDAT)),KWORK(L(LVECDAT)),
     *             KWORK(L(LISTPAR)),DWORK(L(LDSTPAR)),
     *             KWORK(L(LIMGPAR)),DWORK(L(LDMGPAR)),
     *             KWORK(L(LIASMBL)),DWORK(L(LDASMBL)),
     *             LUP,LRHS,NUVP)

C     Now all the memory we can allocate in advance is allocated.
C     Form the handle-list which is passed to the optimizér
C     as integer parameter block.

      HLIST(OLTRIAS ) = LTRIAS
      HLIST(OLCURTR ) = LCURTR
      HLIST(OLMATDAT) = LMATDAT
      HLIST(OLVECDAT) = LVECDAT
      HLIST(OLISTPAR) = LISTPAR
      HLIST(OLDSTPAR) = LDSTPAR
      HLIST(OLIMGPAR) = LIMGPAR
      HLIST(OLDMGPAR) = LDMGPAR
      HLIST(OLIASMBL) = LIASMBL
      HLIST(OLDASMBL) = LDASMBL
      HLIST(OLIGEOM ) = LIGEOM
      HLIST(OLDGEOM ) = LDGEOM
      HLIST(ODUP    ) = LUP
      HLIST(ODRHS   ) = LRHS

      HLIST(ONUVP   ) = NUVP !KWORK(L(LVECDAT)+(NLMAX-1)*SZN2VI+ONEQV-1)
      
      HLIST(OLIOPTP ) = LIOPTP
      HLIST(OLDOPTP ) = LDOPTP

C     Open a file where to write the output to:

      WRITE (CBFN,'(A)') 'optimization.res'
      
      IFMTS = 1
      HLIST(OIFXOUT) = 67
      CALL OF0 (HLIST(OIFXOUT),CBFN,IFMTS)

C     Duplicate the original meshes into CURTRI. Depending on how we
C     deform the mesh, we don't have to duplicate everything...
      
C     For IGAMET<>0 don't duplicate structural information (like KVERT), 
C     as this is not changed in the grid deformation process. Other
C     information like coordinates,... we have to back up so we can
C     restore it when necessary.

      IF (IGAMET.NE.0)
     *  J = 2**31-1 - 2**0 - 2**1 - 2**7 - 2**12 - 2**13 - 2**15

C     For IMESH=0 or if we don't perform optimization, we can share 
C     all information with TRIAS, since nothing will be changed:

      IF ((IGAMET.EQ.0).OR.(IOPTTP.EQ.0)) J = 2**31 - 2**0
     
C     Start the grid duplication for all levels
      
      DO I=NLMAX,NLMIN,-1

C       Copy DCORVG only on the maximum level, since the information
C       is shared between all levels!

        IF (I.EQ.NLMAX) THEN
          CALL TRIDUP(KWORK(L(LTRIAS)+SZTRIA*(I-1)),
     *                KWORK(L(LCURTR)+SZTRIA*(I-1)),0,J)
        ELSE
          CALL TRIDUP(KWORK(L(LTRIAS)+SZTRIA*(I-1)),
     *                KWORK(L(LCURTR)+SZTRIA*(I-1)),0,IOR(J,2**0))
          KWORK(L(LCURTR)+SZTRIA*(I-1)+OLCORVG-1) = 
     *          KWORK(L(LCURTR)+SZTRIA*(NLMAX-1)+OLCORVG-1)
        END IF
        
      END DO

      IF (IOPTPR.EQ.0) THEN

C       Initialize minimum + maximum coordinates, init. step length.
C       The LCOORDS-array in the probing subroutine will have
C       this structure: X-Coordinate, Y-Coordinate, radius.

        DSTART(1) = OBXMIN
        DSTART(2) = OBYMIN
        DSTART(3) = OBRMIN
        
        DSTEP(1) = OBXSTP
        DSTEP(2) = OBYSTP
        DSTEP(3) = OBRSTP

        DEND(1) = OBXMAX
        DEND(2) = OBYMAX
        DEND(3) = OBRMAX
      
      ELSE IF (IOPTPR.EQ.10) THEN
      
C       Initialize min/max coordinates and init. step length.

        NDIM = 2

        DSTART(1) = OBYINI 
C       0.5*(OBYMIN+OBYMAX)
        DSTART(2) = OBYINI+OBYSTP
C        +OBYSTP
C       0.5*(OBYMIN+OBYMAX)
        
        DSTEP(1) = OBYSTP
        DSTEP(2) = OBYSTP

        DEND(1) = OBYMAX
        DEND(2) = OBYMAX
        
      END IF
        
C     Call the optimizer...

      IF (IOPTTP.EQ.0) THEN

C       IOPTTP=0: no optimization, standard CC2D calculation

        CALL CC2DST (MSHOW,HLIST)
        
      ELSE IF (IOPTTP.EQ.1) THEN
      
        PRINT *,'IOPTTP=1 currently not supported.'
        STOP

      ELSE IF (IOPTTP.EQ.2) THEN
      
C       IOPTTP=2: brute force test - all X-,Y-coordinates, 
C                 rotation angles
        
        CALL OPTA00 (NDIM,DSTART,DSTEP,DEND,
     *               HLIST,DLIST,2,DPRSTS,DNXSTS,
     *               DFINAL,STEPS)
      
      ELSE IF (IOPTTP.EQ.100) THEN
      
C       IOPTTP=100: brute force test - as configured above
        
        CALL OPTA00 (NDIM,DSTART,DSTEP,DEND,
     *               HLIST,DLIST,2,DPRSTS,DNXSTS,
     *               DFINAL,STEPS)
      
      ELSE IF (IOPTTP.EQ.101) THEN
      
C       IOPTTP=101: Compass search

        CALL OPTA01 (NDIM,DSTART,DSTEP,DOPEPS,OITMAX,
     *               HLIST,DLIST,2,DPRSTS,DNXSTS,
     *               DFINAL,FINTOL,STEPS)

      ELSE IF (IOPTTP.EQ.102) THEN
      
C       IOPTTP=102: Nelder Mead
        
        CALL OPTA02 (NDIM,DSTART,DSTEP,DOPEPS,OITMAX,
     *               0.5D0,0.5D0,2D0,
     *               HLIST,DLIST,2,DPRSTS,DNXSTS,
     *               DFINAL,FINTOL,STEPS)

      ELSE IF (IOPTTP.EQ.103) THEN
      
C       IOPTTP=103: Gradient search - mini-implementation only.
        
        CALL OPTA03 (1,DSTART(2),DSTEP(2),DOPINS,DOPEPS,OITMAX,
     *               HLIST,DLIST,2,DPRSTS,DNXSTS,
     *               DFINAL,FINTOL,STEPS)

      END IF

C     Close output file for opptimizer results

      CLOSE (HLIST(OIFXOUT))
     
C     Delete copies of the triangulation.
C     Node that this will only delete information that is maintained
C     by the CURTRI structures! Information that is shared between
C     CURTRI and TRIAS will not be deleted!

      DO I=NLMAX,NLMIN,-1
        CALL TRIDEL (KWORK(L(LCURTR)+SZTRIA*(I-1)))
      END DO
      
C     Release matrices
     
      CALL DNSTSL (NLMIN,NLMAX,
     *             KWORK(L(LMATDAT)),KWORK(L(LVECDAT)),
     *             KWORK(L(LISTPAR)),KWORK(L(LDSTPAR)),LUP,LRHS)
     
C     Release the parameter blocks, finish.

      CALL ZDISP(0,LMATDAT, 'MATDAT')
      CALL ZDISP(0,LVECDAT, 'VECDAT')
      CALL ZDISP(0,LISTPAR, 'ISTPAR')
      CALL ZDISP(0,LDSTPAR, 'DSTPAR')
      CALL ZDISP(0,LIMGPAR, 'IMGPAR')
      CALL ZDISP(0,LDMGPAR, 'DMGPAR')
      CALL ZDISP(0,LIASMBL, 'IASMBL')
      CALL ZDISP(0,LDASMBL, 'DASMBL')
      CALL ZDISP(0,LCURTR , 'CURTRI')
      
9000  FORMAT(79('-'))
9001  FORMAT(60('-'))

      END
