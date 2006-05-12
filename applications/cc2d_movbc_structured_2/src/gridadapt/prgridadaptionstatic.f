************************************************************************
* This file contains subroutines to prepare the solver of the
* grid adaption routines of the static method.
*
* All routines base on the parameter blocks TM020IExtParams/
* TM020DExtParams, whose structure is defined in the include file
* "PRGRIDADAPTIONSTATICSOLVER.INC".
*
* The callback routines in this file are specialized to solve a
* pure Neumann problem with the M020 multigrid solver. 
* For this duty they access the information in the TM020IExtParams/
* TM020DExtParams structures, which is passed to them by the solver.
************************************************************************

************************************************************************
************************************************************************

************************************************************************
* Initialise multigrid solver for grid adaption
*
* - Initializes the general multigrid structures
* - Allocates and calculates the structure of the 
*   laplacian=system matrices on all levels,
* - Allocates memory for RHS, defect and solution vectors on all
*   levels.
* - Calculates the number of unknowns on each level.
* - Initialises the handles of the matrix entries and the smoother
*   matrices with 0 to mark them as "not calculated".
* - Calculates the symbolical factorisation of the matrix on the
*   coarse grid if UMFPACK4-solver is used there.
*
* The routine needs a set of grid structures TRIA that contain
* geometry information about all grids that are used in multigrid.
* The array has NNLEV elements, wher only the elements
* NLMIN..NLMAX have to be initialized with data.
*
* In:
*   IPARAM : array [1..SZGLSI] of integer = TIGridAdaptionSolver
*            Integer parameter array to configure the solver
*   DPARAM : array [1..SZGLSD] of double = TDGridAdaptionSolver
*            Double precision parameter array to configure the solver
*   ILEV   : maximum level that should be used in the multigrid.
*   TRIAS  : array [1..NNLEV] of STRIA
*            = array [1..SZTRIA,1..NNLEV] of integer
*            Triangulation structures for every level determining
*            the grids. Needed by the matrix generation routines.
*
* Out:
*   IMGPAR : array [1..SZGSSI] of integer = TM020IExtParams
*            Integer parameter array for multigrid algorithm,
*            filled with data.
*   DMGPAR : array [1..SZGSSD] of double  = TM020DExtParams
*            Double precision parameter array for multigrid algorithm,
*            filled with data.
*   LSOL   : handle to solution vector on finest level
*   LRHS   : handle to RHS-vector on finest level
*   NEQFL  : size of the system on finest level
*
* The structure information of the bilinear forms must be initialised
* with INGSGS before calling this routine!
* The multigrid structure has to be initialized by INM020 before
* calling this routine!
************************************************************************

      SUBROUTINE INMGG2 (IPARAM, DPARAM, ILEV, TRIAS, IMGPAR, DMGPAR,
     *                   LSOL, LRHS, NEQFL)
      
      IMPLICIT NONE
      
      INCLUDE 'cout.inc'
      INCLUDE 'cmem.inc'
      INCLUDE 'cerr.inc'

      INCLUDE 'cbasicelem.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'sgridadaptsolverparams.inc'
      
      INCLUDE 'stria.inc'

      INCLUDE 'ssolvers.inc'
      INCLUDE 'm020.inc'
      INCLUDE 'sgridadaptstaticsolver.inc'
      
      INTEGER IMGPAR(SZGSSI), IPARAM(SZGLSI),ILEV,LSOL,LRHS,NEQFL
      DOUBLE PRECISION DMGPAR(SZGSSD), DPARAM(SZGLSD)
      INTEGER TRIAS(SZTRIA,NNLEV)
      
C externals      
      
      EXTERNAL E011
      
C local variables
      
      INTEGER ISYMM,I,J,LDEGU,LCOLH,NDEGU,STARTI,STARTD
      INTEGER LTMPA,LCLTMP,LLDTMP,NMILV,NMXLV
      INTEGER LA,LCOL,LLD,KTR1,KTR2,NEQ,NA
      INTEGER TRITMP(SZTRIA)
      
C Initialise multigrid parameter arrays with standard values

      CALL LCL3(IMGPAR,SZGSSI)
      CALL LCL1(DMGPAR,SZGSSD)
      CALL INM020(IMGPAR, DMGPAR)
      
C Transfer triangulation structures into multigrid block.
C Necessary for prolongation/restriction, which need information
C about the triangulation.

      CALL LCP3(TRIAS,IMGPAR(OTRIAS),NNLEV*SZTRIA)

C determine minimum and maximum level
      
      NMXLV = ILEV
      
      NMILV = IPARAM(OICGSLV)
      IF (NMILV.LT.1) NMILV=NMXLV+NMILV
      IF (NMILV.LT.1) NMILV=1
      IF (NMILV.GT.NMXLV) NMILV = NMXLV
      
C Set up some parameters of the multigrid solver that are problem-
C specific: We are using filtering as we solve a pure Neumann problem!

      IMGPAR(OFILT) = 1
      
C Transfer the solver parameters from the parameter block IPARAM/DPARAM
C to the solver structure to configure the solver how to react:

      IMGPAR(OSLTAG) = 11
      IMGPAR(ONLMIN) = NMILV
      IMGPAR(ONLMAX) = NMXLV
      IMGPAR(ONITMAX) = IPARAM(ONITSGS)
      
      IF (IPARAM(OIRELGS).EQ.0) THEN
        DMGPAR (OEPSREL) = 0D0
        DMGPAR (OEPSABS) = DPARAM(OEPSGSS)
      ELSE
        DMGPAR (OEPSREL) = DPARAM(OEPSGSS)
        DMGPAR (OEPSABS) = 0D0
      END IF
      
C Now we are going to precalculate as much information as possible,
C more precisely we precalculate everything that does not change during
C the grid adaption process. This includes handles to matrix structure,
C sorting strategy, handles to RHS,...
C All computed information is stored in the extended multigrid
C information block, which is passed to the callback routines during the
C solving process later.
C
C Remember to set the current level back when the subroutine is completed.
C This is important if we are doing multiple grid adaption passes
C on different levels - and anywhere it's cleaner, otherwise we would
C have a side effect...
C So make a backup of the current triangulation.

C      CALL C2TRIA(TRITMP)
      
C Create system/mass matrix structure on all levels
      
      IF (MT.GE.2) THEN
        WRITE(MTERM,'(A$)') 'Creating pointer vectors for '//
     *        'all levels: ['
        !CALL FLUSH (MTERM)
      END IF 
       
C don't exploit any symmetry

      ISYMM=0
      J=1

      DO I=NMILV,NMXLV
        
        IF (MT.GE.2) THEN
          WRITE(MTERM,'(I2$)') I
          !CALL FLUSH(MTERM)
        END IF

C As long as the matrix generation routines rely on parameters in the
C COMMON blocks, we have to copy our grid information from the TRIA-
C structure of the current level to that before calling the generation
C routine. Call TRIA2C to copy the information necessary for matrix
C generation:
C
C        CALL TRIA2C ()
        
C Generate structure of system matrix

        CALL XAP7X(IMGPAR(OLCOL+I-1),IMGPAR(OLLD+I-1),NA,
     *             NEQ,TRIAS(1,I),E011,ISYMM)
        IF (IER.NE.0) GOTO 99998
        
        IMGPAR(OKNA+I-1) = NA
        IMGPAR(OKNEQ+I-1) = NEQ

C Allocate memory for solution/RHS,...
C Write handles into the solver structure.

        CALL ZNEW(NEQ,1,IMGPAR(OLSOL+I-1),'DX    ')
        CALL ZNEW(NEQ,1,IMGPAR(OLRHS+I-1),'DB    ')
        CALL ZNEW(NEQ,1,IMGPAR(OLAUXD+I-1),'DD    ')
        CALL ZNEW(NEQ,1,IMGPAR(OLAUX2+I-1),'DAUX  ')
        IF (IER.NE.0) GOTO 99998

C All other variables initialise with 0 to mark them as "not calculated"

        IMGPAR(OLA+I-1)    = 0
        IMGPAR(OLILU0+I-1) = 0
        IMGPAR(OLTR1+I-1)  = 0
        IMGPAR(OLTR2+I-1)  = 0
        
      END DO

      IF (MT.GE.2) WRITE(MTERM,'(A)') ']'

C Return the handles of the solution vector and the RHS-vector
C on the finest level to the caller, so that the caller can
C write data to and read the solution from it:

      LSOL = IMGPAR(OLSOL+NMXLV-1)
      LRHS = IMGPAR(OLRHS+NMXLV-1)
      NEQFL = IMGPAR(OKNEQ+NMXLV-1)

      IF (IPARAM(OISRTGS).GT.0) THEN
        
C If sorting is active, use the matrix structure to calculate the
C permutations. Don't really resort the matrix yet.

C Cuthill-McKee resorting; use matrix structure to calculate 
C the permutation.
          
        IF (IPARAM(OISRTGS).EQ.1) THEN

          IF (MT.GE.2) THEN
            WRITE(MTERM,'(A$)') 'Calculating renumbering strategy... ['
            !CALL FLUSH (MTERM)
          END IF 

          DO I=NMILV,NMXLV

            IF (MT.GE.2) THEN
              WRITE(MTERM,'(I2$)') I
              !CALL FLUSH (MTERM)
            END IF
          
            LA = IMGPAR(OLA+I-1)
            LCOL = IMGPAR(OLCOL+I-1)
            LLD = IMGPAR(OLLD+I-1)
            NEQ = IMGPAR(OKNEQ+I-1)
            NA = IMGPAR(OKNA+I-1)
          
C Allocate storage for permutation vectors
          
            CALL ZNEW(NEQ, 3,KTR1,'KKTRA1')
            CALL ZNEW(NEQ, 3,KTR2,'KKTRA2')
            IF (IER.NE.0) GOTO 99998

            IMGPAR(OLTR1+I-1) = KTR1
            IMGPAR(OLTR2+I-1) = KTR2
          
C Allocate auxiliary vectors

            CALL ZNEW(NA,-3,LCOLH,'KCOLH ')
            CALL ZNEW(NEQ,-3,LDEGU,'KDEGU ')
            IF (IER.NE.0) GOTO 99998

C Calculate renumbering strategy for system matrix.
C The mass matrix is not resorted, since no system is solved with that!

            NDEGU = IMGPAR(OKNEQ+I-1)
            CALL CUTCE0(KWORK(L(LLD)),KWORK(L(LCOL)),
     *                  KWORK(L(LCOLH)),KWORK(L(LDEGU)),NEQ,NDEGU)

C Calculate permutation vectors KTR1GS and KTR2GS

            CALL CUTCE1(KWORK(L(LLD)),KWORK(L(LCOLH)),NEQ,
     *                  KWORK(L(KTR1)),KWORK(L(KTR2)))

C That's all; free auxiliary vectors

            CALL ZDISP (0,LDEGU,'KDEGU ')
            CALL ZDISP (0,LCOLH,'KCOLH ')
            
          END DO

          IF (MT.GE.2) WRITE(MTERM,'(A)') ']'

C Up to now the matrices are not sorted, but the resorting strategy has
C been set up to Cuthill-McKee; so we set IMGPAR(OISRT) to -1.
C IMGPAR(OISRT) always stores the currently active sorting strategy of the
C matrices. This is necessary because when the matrices are recalculated
C during the iteration, the calculation routines need unsorted matrices.
C So we have to know if we have to sort the matrix structure back
C before recalculation or not.
      
          IMGPAR(OISRT) = -1
            
        END IF
        
      END IF
      
C If we use UMFPACK4 as coarse grid solver, compute the symbolical factorisation

      IF (IPARAM(OICGSGS).EQ.1) THEN
        
        IF (MT.GE.2) THEN
          WRITE(MTERM,'(A$)') 'Calculating sym. factorisation on '//
     *                        'coarse grid... ['
          !CALL FLUSH (MTERM)
        END IF 
        
C First create a copy of our matrix on the coarse grid

        NA = IMGPAR(OKNA+NMILV-1)
        NEQ = IMGPAR(OKNEQ+NMILV-1)
                
        CALL ZNEW (NA,-1,LTMPA,'LTMPA ')
        CALL ZNEW (NA,-3,LCLTMP,'LCLTMP')
        CALL ZNEW (NEQ+1,-3,LLDTMP,'LLDTMP')

C Fill the matrix with 1 to prevent the entries to be treated as 0

        DO I=1,NA
          DWORK(L(LTMPA)+I-1) = 1D0
        END DO

C We have to renumber the copy of our matrix to lower memory usage, 
C and to ensure that in the later solution process the renumbered vector
C fits to the numerically factorised matrix.
C
C If we are not using renumbering, copy the structure.
C If we are using renumbering, use the Cuthill-McKee renumbering routine
C to build the matrix structure:

        KTR1 = IMGPAR(OLTR1+NMILV-1)
        KTR2 = IMGPAR(OLTR2+NMILV-1)
        LCOL = IMGPAR(OLCOL+NMILV-1)
        LLD  = IMGPAR(OLLD+NMILV-1)
      
        IF (IPARAM(OISRTGS).GT.0) THEN
        
C The following call uses source=dest. matrix, which seems suspicious
C but is ok here! LTMP is only a temporary matrix vector for the
C factorisation and is deleted afterwards!
        
          CALL MTSRTD(DWORK(L(LTMPA))  ,DWORK(L(LTMPA)),
     *                KWORK(L(LCLTMP)),KWORK(L(LCOL)),
     *                KWORK(L(LLDTMP)),KWORK(L(LLD)),
     *                KWORK(L(KTR1)),KWORK(L(KTR2)),
     *                NEQ)

        ELSE

          CALL ZCPY (LCOL,'KCOLA ',LCLTMP,'LCLTMP')
          CALL ZCPY (LLD,'KLDA  ',LLDTMP,'LLDTMP')

C the matrix itself is already filled with 1...
        
        END IF

C Resort the entries in the matrix to ascending order for UMF4

        CALL CSRSRT (DWORK(L(LTMPA)),KWORK(L(LCLTMP)),KWORK(L(LLDTMP)),
     *               NEQ)
     
C The column/row indices must be descremented by 1 to be 0-based,
C because UMFPACK4 uses 0-based vectors!

        CALL M7IDSH (KWORK(L(LCLTMP)),KWORK(L(LLDTMP)),NEQ)
     
C Initialize the Umfpack4 control structure 

        CALL UMF4DEF(DMGPAR(OCNTRU4))
        
C Symbolically factorize the matrix; this gives KU4SYM.

        CALL UMF4SYM(NEQ,NEQ,KWORK(L(LLDTMP)),KWORK(L(LCLTMP)),
     *               DWORK(L(LTMPA)),IMGPAR(OKU4SYM),
     *               DMGPAR(OCNTRU4),DMGPAR(OINFOU4))
     
        IF (DMGPAR(OINFOU4).NE.0) THEN
          WRITE (MTERM,'(A,F20.10)') 'Error: ',DMGPAR(OINFOU4)
          IER=1
        END IF
        IF (IER.NE.0) GOTO 99998

C Don't perform a numerical factorisation - this is done later...

C The copy of the matrix is no longer necessary; delete it.
      
        CALL ZDISP (0,LLDTMP,'LLDTMP')
        CALL ZDISP (0,LCLTMP,'LCLTMP')
        CALL ZDISP (0,LTMPA,'LTMPA ')

        IF (MT.GE.2) WRITE(MTERM,'(A)') 'ok]'

      END IF
      
C Coarse grid solver structure setup

      STARTI = OCGSI
      STARTD = OCGSD
      
C The coarse grid solver needs a
C TSolverIMGPARs/TSolverDMGPARs parameter block with information
C how to behave. This parameter block is placed at position
C IMGPAR(DCGSD) / DMGPAR(DCGSD) of the extended multigrid parameter
C block. For the whole grid adaption process there is a
C definite behaviour of that solver, so we can initialize it in
C advance:

      CALL INGSLV (IMGPAR(STARTI),DMGPAR(STARTD))
      
C Transfer the message level of the coarse grid solver from the
C multigrid structure to the coarse-grid-solver structure

C      IMGPAR(OMCGTRM) = 2
      IMGPAR(STARTI+OMSGTRM-1) = IMGPAR(OMCGTRM)
      
C If we are using BiCGStab as coarse-grid solver, we have to
C prepare that!

      IF (IPARAM(OICGSGS).EQ.0) THEN
      
C 1.) The accuracy of the solver has to be configured as defined
C     in the parameters of the DAT file:

        IF (IPARAM(OIRCGGS).NE.0) THEN
          DMGPAR(STARTD+OEPSREL-1) = DPARAM(OEPSCGS)
          DMGPAR(STARTD+OEPSABS-1) = 0D0
        ELSE
          DMGPAR(STARTD+OEPSREL-1) = 0D0
          DMGPAR(STARTD+OEPSABS-1) = DPARAM(OEPSCGS)
        END IF
        IMGPAR(STARTI+ONITMAX-1) = IPARAM(ONICGGS)
      
C 2.) The damping parameter can has to be configured as defined
C     in the parameters of the DAT file. The damping parameter is
C     always saved in the OMEGA-parameter of DMGPAR:

        DMGPAR(STARTD+OOMGPC-1) = DPARAM(OOMCGGS)
      
C 3.) We don't need any timings; this is done by the multigrid
C     algorithm internally

        IMGPAR(STARTI+OTIM-1) = 0
   
C 4.) The coarse grid solver needs filtering to be turned on, since
C     there is a Neumann problem considered on the lowest level:

        IMGPAR(STARTI+OFILT-1) = 1
      
C 5.) We can initialize whether a preconditioner is to be used in the
C     coarse grid solver:

        IF (IPARAM(OIPCGGS).EQ.0) IMGPAR(STARTI+OPCTAG-1) = 1
        IF (IPARAM(OIPCGGS).EQ.1) IMGPAR(STARTI+OPCTAG-1) = 4
        IF (IPARAM(OIPCGGS).EQ.2) IMGPAR(STARTI+OPCTAG-1) = 8
      
C     ... but if ILU is used as preconditioner, we can't initialize this
C     here, as the handle is still missing. So this is postponed until
C     we have that handle.

C Write into that solver structure the type of the coarse grid solver
C we are using. We use the solver tag to save that information.

        IMGPAR(STARTI+OSLTAG-1) = 6

      ELSE

C Write into that solver structure that we use UMFPACK4 as coarse
C grid solver:

        IMGPAR(STARTI+OSLTAG-1) = 13
        
C For the UMFPACK4 solver we don't use an extra solver block as
C it's done for the BiCGStab-solver above. This type of solver
C is directly handled by the callback-routine. All necessary
C information for it is calculated later in the preparation routine
C below and saved in a special block in the double-precision 
C parameter array.

      END IF
      
C Smoother solver structure setup

      STARTI = OSMSI
      STARTD = OSMSD

C The smoother needs a TSolverIMGPARs/TSolverDMGPARs 
C parameter block with information how to behave. 
C This parameter block is placed at position
C IMGPAR(DSMSD) / DMGPAR(DSMSD) of the extended multigrid parameter
C block. For the whole grid adaption process there is a
C definite behaviour of that solver, so we can initialize it in
C advance:

      CALL INGSLV (IMGPAR(STARTI),DMGPAR(STARTD))

C Mark which type of smoother we are using - in the SLTAG variable
C of the smoother structure. The callback-routine will use this value
C later to decide how to react. Furthermore define the SMTAG-variable
C in the multigrid structure to indicate the type of the smoother.
C
C For the moment we only have two possibilities: no smoother (PCSMG=-1)
C and BiCGStab-smoother (PCSMG>=0):

      IF (IPARAM(OIPCSMG).EQ.0) THEN
        IMGPAR(STARTI+OSLTAG-1) = 206
        IMGPAR(OSMTAG) = 6
      END IF

C Transfer the message level of the smoother from the
C multigrid structure to the smoother structure

      IMGPAR(STARTI+OMSGTRM-1) = IMGPAR(OMSMTRM)
      
C In the other case both variables stay 0 to mark that no smoother
C is used at all.

C The damping parameter can has to be configured as defined
C in the parameters of the DAT file. We can save it either in the
C smoother structure as OMEGA, or in the OMGSM-variable of the
C multigrid structure. We only have to make sure that the callback
C routine finds the value again. As we use one parameter on all levels
C and as the smoother structure is more or less temporary (some
C variables are changed in the callback routine to fit the current
C level), we choose to save the damping parameter in the OMGSM-
C variable of DMGPAR:

      DMGPAR(OOMGSM) = DPARAM(OOMGPGS)
      
C Also the smoother should use filtering:

      IMGPAR(STARTI+OFILT-1) = 1
      
C If we are using BiCGStab as smoother, we have to prepare more!

      IF (IPARAM(OIPCSMG).NE.-1) THEN
      
C     We can initialize whether a preconditioner is to be used in the
C     BiCGStab-smoother. The IPCSGS-parameter defines the preconditioner

        IF (IPARAM(OIPCSGS).EQ.0) IMGPAR(STARTI+OPCTAG-1) = 1
        IF (IPARAM(OIPCSGS).EQ.1) IMGPAR(STARTI+OPCTAG-1) = 4
        IF (IPARAM(OIPCSGS).EQ.2) IMGPAR(STARTI+OPCTAG-1) = 8
      
C     ... but if ILU is used as preconditioner for the smoother, we 
C     can't initialize this here, as the handles to the matrix are still
C     missing. So this is postponed until we have that handle.

      END IF

C Restore the original triangulation

99998 CONTINUE
C     CALL TRIA2C(TRITMP)
      
      END

************************************************************************
* Prepare multigrid solver for grid adaption
*
* The structure of the systems have to be initialised with INMGGS
* before!
*
* -Calculates the system matrices on all levels
* -Calculates the ILU0-smoother on all levels if necessary
* -Calculates the reference right-hand-side LONEGS
* -Initializes the KOFFX/KOFFB/KOFFD index arrays for the solver,
*  initializes the number of smoothing steps on each level
* 
* If the matrices/RHS already exist, they will be overwritten.
*
* The routine needs a set of grid structures TRIA that contain
* geometry information about all grids that are used in multigrid.
* The array has NNLEV elements, wher only the elements
* NLMIN..NLMAX have to be initialized with data.
*
* After the routine is finished, all data for the multigrid solver
* is initialized. M020 can - and should - be called directly. 
* To prevent shifts in the L-array, which would destroy the
* correct initialization of the multigrid solver, no arrays
* must be deleted between this function and the return of
* the M020 solver!
*
* In:
*   IPARAM : array [1..SZGLSI] of integer = TIGridAdaptionSolver
*            Integer parameter array to configure the solver
*   DPARAM : array [1..SZGLSD] of double = TDGridAdaptionSolver
*            Double precision parameter array to configure the solver
*   IMGPAR : array [1..SZGSSI] of integer
*            Integer parameter array for multigrid algorithm,
*            filled with data, prepared by the initialization routine.
*   DMGPAR : array [1..SZGSSD] of double
*            Double precision parameter array for multigrid algorithm,
*            filled with data, prepared by the initialization routine.
*   TRIAS  : array [1..NNLEV] of STRIA
*            = array [1..SZTRIA,1..NNLEV] of integer
*            Triangulation structures for every level determining
*            the grids. Needed by the matrix generation routines.
*   ILEV   : Maximum level that should be used in multigrid
* Out:
*   The handles in IMGPAR/DMGPAR are filled with data.
*   Smoother- and coarse-grid solver structures are prepared
*   for being used.
************************************************************************

      SUBROUTINE PRMGG2 (IPARAM, DPARAM, TRIAS, IMGPAR, DMGPAR)
      
      IMPLICIT NONE
      
      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      INCLUDE 'cmem.inc'

      INCLUDE 'cbasicelem.inc'

      INCLUDE 'cbasicmg.inc'
      INCLUDE 'sgridadaptsolverparams.inc'
      
      INCLUDE 'stria.inc'

      INCLUDE 'ssolvers.inc'
      INCLUDE 'm020.inc'
      INCLUDE 'sgridadaptstaticsolver.inc'
      
      INTEGER IMGPAR(SZGSSI),IPARAM(SZGLSI)
      DOUBLE PRECISION DMGPAR(SZGSSD),DPARAM(SZGLSD)
      INTEGER TRIAS(SZTRIA,NNLEV)
      
C externals
      
      EXTERNAL E011,COEAGS,COEMGS
      
C Names of matrix and vector (for messages only)
      
      CHARACTER ARRDA*6
      DATA ARRDA/'DA    '/
      SAVE ARRDA

      CHARACTER ARRDM*6
      DATA ARRDM/'DM    '/
      SAVE ARRDM

C local variables

      INTEGER ISYMM,I,ICLR,LCOLH,LLDH,LAH,STARTI,STARTD
      DOUBLE PRECISION TOLILU,ALPILU
      INTEGER INDILU
      INTEGER LTMPA,LCLTMP,LLDTMP,NMXLV,NMILV
      INTEGER NA,NEQ,LA,LCOL,LLD,KTR1,KTR2,LAILU
      INTEGER TRITMP(SZTRIA)

C Variables to hold the structure of the bilinear forms:

      INTEGER NBLOC
      PARAMETER (NBLOC=1)

      INTEGER KABA(2,NNAB,NBLOC),KABM(2,NNAB,NBLOC),KF(NNAB,NBLOC)
      INTEGER KABAN(NBLOC),KABMN(NBLOC),KFN(NBLOC)
      LOGICAL BCONA(NBLOC),BCONM(NBLOC),BCONF(NBLOC)
      LOGICAL BSNGLA(NBLOC),BSNGLM(NBLOC),BSNGLF(NBLOC)
      
C Handle all levels - up to the current one!
C Remember to set the current level back when the subroutine is completed.
C This is important if we are doing multiple grid adaption passes
C on different levels - and anywhere it's cleaner, otherwise we would
C have a side effect...

      NMXLV = IMGPAR(ONLMAX)
      NMILV = IMGPAR(ONLMIN)
      
C Make a backup of the current triangulation structure

C      CALL C2TRIA(TRITMP)

C If the matrices are to be recalculated, they are perhaps sorted.
C In this case we have to sort the matrix structure back before
C recalculation. After the calculation is completed, we can again
C apply the sorting.
C We only reconstruct the old matrix structure without sorting
C the entries back, as the entries are recalculated=overwritten
C below anyway.

      IF (IMGPAR(OISRT).GT.0) THEN
        IF (MT.GE.2) THEN
          WRITE(MTERM,'(A$)') 'Sorting back matrix structures ['
          !CALL FLUSH (MTERM)
        END IF
        
        DO I=NMILV,NMXLV
        
          IF (MT.GE.2) THEN
            WRITE(MTERM,'(I2$)') I
            !CALL FLUSH(MTERM)
          END IF
      
C Allocate auxiliary vectors that take a copy of the structure

          NA   = IMGPAR(OKNA+I-1)
          NEQ  = IMGPAR(OKNEQ+I-1)
          LCOL = IMGPAR(OLCOL+I-1)
          LLD  = IMGPAR(OLLD+I-1) 
          KTR1 = IMGPAR(OLTR1+I-1) 
          KTR2 = IMGPAR(OLTR2+I-1) 

          CALL ZNEW(NA  ,-3,LCOLH,'KCOLH ')
          CALL ZNEW(NEQ+1,-3,LLDH ,'KLDH  ')
          IF (IER.NE.0) GOTO 99998

          CALL ZCPY(LCOL,'KCOLA ',LCOLH,'KCOLH ')
          CALL ZCPY(LLD ,'KLDA  ',LLDH, 'KLDH  ')

C Call the routine to sort the matrix back.
C Take the aux. vectors as source, overwrite the structure
C of the real matrix.

          CALL MTSRTR (KWORK(L(LCOL)),KWORK(L(LCOLH)),
     *                 KWORK(L(LLD)),KWORK(L(LLDH)),
     *                 KWORK(L(KTR1)),KWORK(L(KTR2)),
     *                 NEQ)
     
C Release aux. vectors

          CALL ZDISP (0,LLDH ,'KLDH  ')
          CALL ZDISP (0,LCOLH,'KCOLH ')

        END DO

        IF (MT.GE.2) WRITE(MTERM,'(A)') ']'
        
C Matrix structures are now unsorted, but sorting strategy is
C still prepared. Switch the sign of the sorting-strategy tag 
C to mark the matrices as unsorted.
        
        IMGPAR(OISRT) = -IMGPAR(OISRT)
      END IF

      IF (MT.GE.2) THEN
        WRITE(MTERM,'(A$)') 'Calculating matrices on all '//
     *                      'levels: ['
      END IF

C Initialize bilinear forms.

      CALL INDCGS (7, BCONA, BCONM, BCONF,
     *             BSNGLA, BSNGLM, BSNGLF,
     *             KABAN, KABMN, KFN,
     *             KABA, KABM, KF)

C Don't exploit any symmetry, delete old matrix if there is one

      ISYMM = 0
      ICLR  = 1

      DO I=NMILV,NMXLV
        
C As long as the matrix generation routines rely on parameters in the
C COMMON blocks, we have to copy our grid information from the TRIA-
C structure of the current level to that before calling the generation
C routine. Call TRIA2C to copy the information necessary for matrix
C generation:
C
C        CALL TRIA2C (TRIAS(1,I))
        
        IF (MT.GE.2) THEN
          WRITE(MTERM,'(I2$)') I
          !CALL FLUSH(MTERM)
        END IF
        
        NA   = IMGPAR(OKNA+I-1)
        NEQ  = IMGPAR(OKNEQ+I-1)
        LCOL = IMGPAR(OLCOL+I-1)
        LLD  = IMGPAR(OLLD+I-1) 

C Generate system matrix; resulting handle for the new matrix
C is saved directly into the MG-structure.
        
        CALL XAB7X(IMGPAR(OLA+I-1),LCOL,LLD,NA,NEQ,
     *             NBLOC,ICLR,TRIAS(1,I),E011,.FALSE.,COEAGS,
     *             BCONA,KABA,KABAN,IPARAM(OICUBGS),ISYMM,ARRDA,
     *             0,0D0)

        IF (IER.NE.0) GOTO 99998
      
      END DO

      IF (MT.GE.2) WRITE(MTERM,'(A)') ']'

C Apply the (precalculated) resorting strategy for all system matrices:

      IF (IPARAM(OISRTGS).GT.0) THEN
      
        IF (MT.GE.2) THEN
          WRITE(MTERM,'(A$)') 'Sorting all levels: ['
          !CALL FLUSH (MTERM)
        END IF

        DO I=NMILV,NMXLV
        
          IF (MT.GE.2) THEN
            WRITE(MTERM,'(I2$)') I
            !CALL FLUSH(MTERM)
          END IF
      
          NA   = IMGPAR(OKNA+I-1)
          NEQ  = IMGPAR(OKNEQ+I-1)
          LA   = IMGPAR(OLA+I-1)
          LCOL = IMGPAR(OLCOL+I-1)
          LLD  = IMGPAR(OLLD+I-1) 
          KTR1 = IMGPAR(OLTR1+I-1) 
          KTR2 = IMGPAR(OLTR2+I-1) 
      
C Allocate auxiliary vectors that take a copy of the matrix

          CALL ZNEW(NA  ,-1,LAH  ,'DAH   ')
          CALL ZNEW(NA  ,-3,LCOLH,'KCOLH ')
          CALL ZNEW(NEQ+1,-3,LLDH ,'KLDH  ')
          IF (IER.NE.0) GOTO 99998

          CALL ZCPY(LA  ,'DA    ',LAH  ,'DAH   ')
          CALL ZCPY(LCOL,'KCOLA ',LCOLH,'KCOLH ')
          CALL ZCPY(LLD ,'KLDA  ',LLDH, 'KLDH  ')

C Call the routine to sort the matrix.
C Take the aux. vectors as source, overwrite the real matrix.

          CALL MTSRTD(DWORK(L(LA))  ,DWORK(L(LAH)),
     *                KWORK(L(LCOL)),KWORK(L(LCOLH)),
     *                KWORK(L(LLD)),KWORK(L(LLDH)),
     *                KWORK(L(KTR1)),KWORK(L(KTR2)),
     *                NEQ)
     
C Release aux. vectors

          CALL ZDISP (0,LLDH ,'KLDH  ')
          CALL ZDISP (0,LCOLH,'KCOLH ')
          CALL ZDISP (0,LAH  ,'KAH   ')
        END DO

        IF (MT.GE.2) WRITE(MTERM,'(A)') ']'
        
C Matrices are sorted, remember this fact for next recalculation.
C Simply switch the sign again to indicate that the matrices are again
C sorted.
        
        IMGPAR(OISRT) = -IMGPAR(OISRT)
      END IF

C Calculate ILU(0)-matrix where necessary...
     
      IF ((IPARAM(OIPCSGS).EQ.2).OR.(IPARAM(OIPCGGS).EQ.2)) THEN
      
        IF (MT.GE.2) THEN
          WRITE(MTERM,'(A$)') 'Calculating ILU0-matrices... ['
          !CALL FLUSH (MTERM)
        END IF

        TOLILU=1D-15
        ALPILU=0.0D0
        INDILU=1

C On the coarse grid compute the ILU-matrix only if necessary for 
C the solver
        
        IF (IPARAM(OIPCGGS).EQ.2) THEN

          IF (MT.GE.2) THEN
            WRITE(MTERM,'(I2$)') NMILV
            !CALL FLUSH (MTERM)
          END IF

          NA   = IMGPAR(OKNA+NMILV-1)
          NEQ  = IMGPAR(OKNEQ+NMILV-1)
          LA   = IMGPAR(OLA+NMILV-1)
          LCOL = IMGPAR(OLCOL+NMILV-1)
          LLD  = IMGPAR(OLLD+NMILV-1) 

C Duplicate the system matrix

C Allocate a new matrix, if we don't have one yet

          IF (IMGPAR(OLILU0+NMILV-1).EQ.0) THEN
            CALL ZNEW(NA,-1,IMGPAR(OLILU0+NMILV-1),'DAILU ')
            IF (IER.NE.0) GOTO 99998
          END IF
          
          LAILU = IMGPAR(OLILU0+NMILV-1) 
          
C Copy our current matrix to DAILU to serve as source matrix
          
          CALL ZCPY(LA,'DA    ',LAILU,'DAILU ')

C In-Place replace it by the ILU(0)-matrix

          CALL IFD17(DWORK(L(LAILU)),KWORK(L(LCOL)),
     *             KWORK(L(LLD)),NEQ,
     *             INDILU,ALPILU,TOLILU)
     
CDEBUG
C            FMT (1) = '(D25.10)'
C            CALL XOWM7 (KAILU0(NMILV),KCLAGS(NMILV),KLDAGS(NMILV),
C     *                  KNEQGS(NMILV),
C     *                  'MATRIX',70,'#gmv/mat.txt',.true.)
C            CLOSE (70)
CDEBUG

          IF (IER.NE.0) GOTO 99998
          
        END IF

C On the other grids compute it for the smoother

        IF (IPARAM(OIPCSGS).GT.0) THEN
        
          DO I=NMILV+1,NMXLV
          
            IF (MT.GE.2) THEN
              WRITE(MTERM,'(I2$)') I
              !CALL FLUSH (MTERM)
            END IF
            
            NA   = IMGPAR(OKNA+I-1)
            NEQ  = IMGPAR(OKNEQ+I-1)
            LA   = IMGPAR(OLA+I-1)
            LCOL = IMGPAR(OLCOL+I-1)
            LLD  = IMGPAR(OLLD+I-1) 
            KTR1 = IMGPAR(OLTR1+I-1) 
            KTR2 = IMGPAR(OLTR2+I-1) 
            
C Duplicate the system matrix

C Allocate a new matrix, if we don't have one yet

            IF (IMGPAR(OLILU0+I-1).EQ.0) THEN
              CALL ZNEW(NA,-1,IMGPAR(OLILU0+I-1),'DAILU ')
              IF (IER.NE.0) GOTO 99998
            END IF
  
            LAILU = IMGPAR(OLILU0+I-1)

C Copy our current matrix to DAILU to serve as source matrix
            
            CALL ZCPY(LA,'DA    ',LAILU,'DAILU ')

C In-Place replace it by the ILU(0)-matrix

            CALL IFD17(DWORK(L(LAILU)),KWORK(L(LCOL)),
     *           KWORK(L(LLD)),NEQ,INDILU,ALPILU,TOLILU)

            IF (IER.NE.0) GOTO 99998
            
          END DO
        END IF
        
        IF (MT.GE.2) WRITE(MTERM,'(A)') ']'

      END IF

C If we use UMFPACK4 as coarse grid solver, compute the numerical
C factorisation; the symbolical factorisation already exists - this
C was computed in the initialization routine.

      IF (IPARAM(OICGSGS).EQ.1) THEN
        
        IF (MT.GE.2) THEN
          WRITE(MTERM,'(A$)') 'Calculating num. factorisation on '//
     *                        'coarse grid... ['
          !CALL FLUSH (MTERM)
        END IF 

C First create a copy of our matrix on the coarse grid

        NA   = IMGPAR(OKNA+NMILV-1)
        NEQ  = IMGPAR(OKNEQ+NMILV-1)
        LA   = IMGPAR(OLA+NMILV-1)
        LCOL = IMGPAR(OLCOL+NMILV-1)
        LLD  = IMGPAR(OLLD+NMILV-1) 
                
        CALL ZNEW (NA,-1,LTMPA,'LTMPA ')
        CALL ZNEW (NA,-3,LCLTMP,'LCLTMP')
        CALL ZNEW (NEQ+1,-3,LLDTMP,'LLDTMP')
        
C Since the matrix structure is now sorted (in contrast to the
C initialisation routine), we can simply copy
C the structure to our temporary vectors

        CALL ZCPY (LA,'DA    ',LTMPA,'LTMPA ')
        CALL ZCPY (LCOL,'KCOLA ',LCLTMP,'LCLTMP')
        CALL ZCPY (LLD,'KLDA  ',LLDTMP,'LLDTMP')

C Resort the entries in the matrix to ascending order for UMF4

        CALL CSRSRT (DWORK(L(LTMPA)),KWORK(L(LCLTMP)),KWORK(L(LLDTMP)),
     *               NEQ)
     
C The column/row indices must be descremented by 1 to be 0-based,
C because UMFPACK4 uses 0-based vectors!

        CALL M7IDSH (KWORK(L(LCLTMP)),KWORK(L(LLDTMP)),NEQ)
     
C Numerically factorize; this gives KU4NUM. Use the symbolical
C factorisation object from the initialisation.
        
        CALL UMF4NUM(KWORK(L(LLDTMP)),KWORK(L(LCLTMP)),DWORK(L(LTMPA)),
     *               IMGPAR(OKU4SYM),IMGPAR(OKU4NUM),
     *               DMGPAR(OCNTRU4),DMGPAR(OINFOU4))
        
        IF (DMGPAR(OINFOU4).EQ.1) THEN
        
C BAAAAAD, UMFPACK can't factorise our - nearly - singular matrix :(
C (Remember: we have pure Neumann boundary...)
C Small trick: add 1D0 to the first element of the matrix to make it
C nonsingular.
C THIS TRICK SHOULD ONLY BE DONE IN THE COARSE-GRID SOLVER IN A REAL
C MULTIGRID SOLVING PROCESS, otherwise it would make the solution wrong!
        
          DWORK(L(LTMPA))=DWORK(L(LTMPA))+1D0
          CALL UMF4NUM(KWORK(L(LLDTMP)),KWORK(L(LCLTMP)),
     *                 DWORK(L(LTMPA)),IMGPAR(OKU4SYM),IMGPAR(OKU4NUM),
     *                 DMGPAR(OCNTRU4),DMGPAR(OINFOU4))

        END IF

        IF (DMGPAR(OINFOU4).NE.0) THEN
          WRITE (MTERM,'(A,F20.10)') 'Error: ',DMGPAR(OINFOU4)
          IER=1
          STOP
        END IF

C The copy of the matrix is no longer necessary; delete it.
      
        CALL ZDISP (0,LLDTMP,'LLDTMP')
        CALL ZDISP (0,LCLTMP,'LCLTMP')
        CALL ZDISP (0,LTMPA,'LTMPA ')
      
        IF (MT.GE.2) WRITE(MTERM,'(A)') 'ok]'
      END IF

C Now we have to do some preparations for the coarse-grid solver.

      STARTI = OCGSI
      STARTD = OCGSD
      
C The basic coarse-grid-solver structure was initialized 
C in INMGG2, but the matrix handles have not been initialized 
C yet (as he haven't had that). Now we have the handle - 
C for the system matrix as well as (if necessary) for the ILU-matrix.
C
C So now it's the time to initialize the preconditioner in advance,
C if one is used. As Jacobi, SSOR and ILU0 are standard 
C preconditioners, they expect handles of the matrix at index
C PRECI/PRECD of "their" solver structure. The handles have been
C calculated above - we can simply copy them into where the
C callback routines expect them.
C
C At first initialize the matrix handle structure with the handle of
C the matrix on the coarse grid

      LA   = IMGPAR(OLA+NMILV-1)
      LCOL = IMGPAR(OLCOL+NMILV-1)
      LLD  = IMGPAR(OLLD+NMILV-1) 

      IMGPAR(STARTI+OMATI-1)   = LA
      IMGPAR(STARTI+OMATI+1-1) = LCOL
      IMGPAR(STARTI+OMATI+2-1) = LLD
    
C Then initialize the preconditioner, if necessary.
C
C Remember that the the coarse grid solver is called with 
C IPARAM/DPARAM pointing to the solver structure of the coarse 
C grid solver, which points to IMGPAR(OCGSI) /
C DMGPAR(OCGSD)! (as a coarse grid solver is a standard solver...)
C Therefore the callback routine of the coarse grid solver (especiallly
C the preconditioner of the coarse grid solver) only have
C access to the solver structure of the coarse grid solver, not
C to all the multigrid information!
C For that reason all multigrid data necessary in the callback-
C routines of the coarse grid solver has to be noted in the
C coarse-grid solver structure, too!

      IF (IPARAM(OICGSGS).EQ.0) THEN

C For Jacobi and SSOR the preconditioner needs the matrix on the 
C coarsest level; these have been found out some lines above.
C
C In contrast an ILU0-preconditioner needs the ILU-matrix on the 
C coarsest level:

        IF (IPARAM(OIPCGGS).EQ.2) THEN
          LA   = IMGPAR(OLILU0+NMILV-1)
          LCOL = IMGPAR(OLCOL+NMILV-1)
          LLD  = IMGPAR(OLLD+NMILV-1) 
        END IF

        IMGPAR(STARTI+OPRECI-1) = LA
        IMGPAR(STARTI+OPRECI+1-1) = LCOL
        IMGPAR(STARTI+OPRECI+2-1) = LLD
        
      END IF

C A similar initialization has to be performed, if we are using 
C BiCGStab as a smoother; we have to prepare the smoother-structure
C for that. But this is not done here in advance here, but in the
C callback-routine for the smoother. This routine copies the
C handles of the matrix that is just needed before calling the
C actual BiCGStab smoother.

C Prepare the solver to solve the system.
C Set up index vectors and number of smoothing steps on all levels.

      CALL PRM020 (IMGPAR, DMGPAR,NMILV,NMXLV, 
     *             IPARAM(ONSMSGS), IPARAM(ONSMSGS),
     *             IMGPAR(OLSOL), IMGPAR(OLRHS), IMGPAR(OLAUXD))
        
C Finally restore the original triangulation

99998 CONTINUE
C     CALL TRIA2C (TRITMP)

C Since we use pure Neumann-boundary-conditions problem here, 
C we don't have to incorporate any boundary conditions into the 
C matrices, the RHS or anywhere else!

      END

************************************************************************
* Dispose the memory that was used for solving the linear systems
* in the grid adaption with the multigrid solver
*
* - Deletes all vectors that have been allocated in INMGGS and PRMGGS
*
* In:
*   ILEV   : Maximum level that was used in multigrid
*   IPARAM : array [1..SZGLSI] of integer = TIGridAdaptionSolver
*            Integer parameter array to configure the solver
*   DPARAM : array [1..SZGLSD] of double = TDGridAdaptionSolver
*            Double precision parameter array to configure the solver
*   IMGPAR : array [1..SZGSSI] of integer
*            Integer parameter array for multigrid algorithm,
*            filled with data, prepared by the initialization and
*            preparation routine.
*   DMGPAR : array [1..SZGSSD] of double
*            Double precision parameter array for multigrid algorithm,
*            filled with data, prepared by the initialization and
*            preparation routine.
* 
* Out:
*   The dynamic memory that was used by the IMGPAR/DMGPAR structures
*   was released. The appropriate handles in IMGPAR/DMGPAR are
*   set to 0.
************************************************************************

      SUBROUTINE DISGS2(IPARAM,DPARAM, IMGPAR,DMGPAR)
      
      IMPLICIT NONE
      
      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      INCLUDE 'cmem.inc'

      INCLUDE 'cbasicelem.inc'

      INCLUDE 'cbasicmg.inc'
      INCLUDE 'sgridadaptsolverparams.inc'
      
      INCLUDE 'stria.inc'

      INCLUDE 'ssolvers.inc'
      INCLUDE 'm020.inc'
      INCLUDE 'sgridadaptstaticsolver.inc'

      INTEGER IMGPAR(SZGSSI),IPARAM(SZGLSI)
      DOUBLE PRECISION DMGPAR(SZGSSD),DPARAM(SZGLSD)

C local variables

      INTEGER I,NMXLV,NMILV
      INTEGER KTR1, KTR2, LAILU, LA

C Handle all levels - up to the current one!
C Obtain the level range from the multigrid parameter block:

      NMXLV = IMGPAR(ONLMAX)
      NMILV = IMGPAR(ONLMIN)

C On all levels delete the memory used. For speed reason do this
C from last level to first level

      DO I=NMXLV,NMILV,-1
      
C some variables might be not initialised:

        KTR1 = IMGPAR(OLTR1+I-1) 
        KTR2 = IMGPAR(OLTR2+I-1) 
        LAILU = IMGPAR(OLILU0+I-1) 
        LA = IMGPAR(OLA+I-1) 

        IF (KTR2.NE.0) CALL ZDISP (0,KTR2,'KTRA2  ')
        IF (KTR1.NE.0) CALL ZDISP (0,KTR1,'KTRA1  ')
        IF (LAILU.NE.0) CALL ZDISP (0,LAILU,'KAILU0 ')
        IF (LA.NE.0) CALL ZDISP (0,LA,'DA    ')
        
        CALL ZDISP (0,IMGPAR(OLAUX2+I-1),'DAUX  ')
        CALL ZDISP (0,IMGPAR(OLAUXD+I-1),'DD    ')
        CALL ZDISP (0,IMGPAR(OLRHS+I-1),'DB    ')
        CALL ZDISP (0,IMGPAR(OLSOL+I-1),'DX    ')
        
        CALL ZDISP (0,IMGPAR(OLLD+I-1),'KLDA  ')
        CALL ZDISP (0,IMGPAR(OLCOL+I-1),'KCLAGS  ')
      END DO

C If necessary, remove UMFPACK4-objects

      IF (IPARAM(OICGSGS).EQ.1) THEN
        IF (IMGPAR(OKU4NUM).NE.0) CALL UMF4FNUM (IMGPAR(OKU4NUM))
        CALL UMF4FSYM (IMGPAR(OKU4SYM))
      END IF

      END
      
************************************************************************
* Prepare right hand side from piecewise linear function
*
* This routine transforms a piecewise linear function into
* a right-hand-side vector of the weak formulation of a Laplacian
* problem -Laplace(u) = f.
*
* KAUX is the starting address of the source function in DWORK, given as
* node values in the NVT nodes of the triangulation TRIA. 
* KRHS is the starting address of the destination vector in DWORK that
* receives the computed right hand side.
*
* PRRLG2 calculates the RHS approximatively using a multiplication
* of a Finite-Difference-like RHS-vector with a mass matrix. The mass
* matrix is build using the cubature formula defined by abs(ICURGS) in
* the parameter block. 
*
* In:
*   IPARAM : array [1..SZGLSI] of integer = TIGridAdaptionSolver
*            Integer parameter array to configure the solver
*   DPARAM : array [1..SZGLSD] of double = TDGridAdaptionSolver
*            Double precision parameter array to configure the solver
*   IMGPAR : array [1..SZGSSI] of integer
*            Integer parameter array for multigrid algorithm,
*            filled with data, prepared by the initialization and
*            preparation routine.
*   DMGPAR : array [1..SZGSSD] of double
*            Double precision parameter array for multigrid algorithm,
*            filled with data, prepared by the initialization and
*            preparation routine.
*   TRIA   : array [1..SZTRIA] of integer
*            Triangulation information of the current grid
*   LB     : Handle to a RHS-vector. The RHS-vector must be
*            given as a vector of function values in the vertices
*            of the triangulation. It's interpreted with the E011
*            element as a piecewise linear function.
*   LRHS   : Handle to a RHS-vector that receives the computed
*            right hand side
*   NEQ    : length of vector identified by LRHS
* Out:
*   The vector specified by LRHS is filled with data.
************************************************************************

      SUBROUTINE PRRLG2 (IPARAM,DPARAM,TRIA,LB,LRHS,NEQ)

      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'cerr.inc'
      INCLUDE 'stria.inc'
      INCLUDE 'cbasicelem.inc'

      INCLUDE 'cbasicmg.inc'
      INCLUDE 'ssolvers.inc'
      INCLUDE 'm020.inc'
      INCLUDE 'sgridadaptstaticsolver.inc'

      INCLUDE 'sgridadaptsolverparams.inc'
      
      INTEGER IPARAM(SZGLSI)
      DOUBLE PRECISION DPARAM(SZGLSD)

      INTEGER NEQ, LB, LRHS, TRIA(SZTRIA)

C Variables to hold the structure of the bilinear forms:

      INTEGER NBLOC
      PARAMETER (NBLOC=1)

      INTEGER KABA(2,NNAB,NBLOC),KABM(2,NNAB,NBLOC),KF(NNAB,NBLOC)
      INTEGER KABAN(NBLOC),KABMN(NBLOC),KFN(NBLOC)
      LOGICAL BCONA(NBLOC),BCONM(NBLOC),BCONF(NBLOC)
      LOGICAL BSNGLA(NBLOC),BSNGLM(NBLOC),BSNGLF(NBLOC)

C externals

      EXTERNAL E011,COEMGS
      
C local variables

      INTEGER ISYMM,ICLR,KCOLD,NVT,I,LCOLD,IC,KRHS,KB,TRIOLD(SZTRIA)
      INTEGER LM,LCOL,LLD,NA,NEQM
      INTEGER TRI2(SZTRIA+1)
      
C Names of matrix and vector (for messages only)
      
      CHARACTER ARRDM*6
      DATA ARRDM/'DM    '/
      SAVE ARRDM

C Copy data from the TRIA-structure for easier access

      NVT = TRIA(ONVT)

C Initialize bilinear forms.

      CALL INDCGS (7, BCONA, BCONM, BCONF,
     *             BSNGLA, BSNGLM, BSNGLF,
     *             KABAN, KABMN, KFN,
     *             KABA, KABM, KF)

C Make a backup of the current triangulation and initialize with
C our current one:

C      CALL C2TRIA(TRIOLD)
C      CALL TRIA2C(TRIA)

C Our right hand side is given as values in the node of the grid.
C To create a FEM-RHS from that, we have to multiply with the mass
C matrix, that is computed if we treat the given vector as piecewise
C linear interpolant of node values.
C
C So we have to create the mass matrix and multiply by it.
C
C For testing purposes we have two cases here, dependent on ICURGS:
C If ICURGS=2, we create a lumped mass matrix and use that one.
C Otherwise we create a real mass matrix with ICURHS as cubature 
C formula and multiply by that one.

      IF (ABS(IPARAM(OICURGS)).EQ.2) THEN

C The easiest mass matrix can be obtained by generating a lumped
C mass matrix. This has exactly NVT elements, which can be stored
C temporarily in the target vector identified by LRHS. Multiplying
C this vector componentwise with the vector identified by KB generates
C the actual RHS vector.
C
C We want to use the XAB07 routine for building this matrix.
C This routine needs a KCOL and KLD vector. As we know our matrix
C is lumped, we can directly set up these vectors, as they are only
C an enumeration 1,2,3,...!
C We reserve one vector for both, KCOL and KLD, because in this special
C case they are the same. Don't forget that KLD has one element more!

        CALL ZNEW (NEQ+1,-3,LCOLD,'KCOLD ')
        KCOLD = L(LCOLD)
        DO I=0,NEQ
          KWORK(KCOLD+I) = I+1
        END DO
        
C Ok, call XAB07 to calculate the matrix. Use KCOL=KLD because in
C this special case they are the same! Use trapezoidal rule for
C calculating the entries to make sure the matrix is lumped.

        IC = 2
        ISYMM = 0
        ICLR = 1

        CALL XAB7X(LRHS,LCOLD,LCOLD,NVT,NVT,
     *             NBLOC,ICLR,TRIA,E011,.FALSE.,COEMGS,
     *             BCONM,KABM,KABMN,IC,ISYMM,ARRDM,0,0D0)

C Restore the triangulation
C
C        CALL TRIA2C (TRIOLD)

C Now multiply that mass matrix with the right hand side vector
C element by element, this gives the actual right hand side:

        KRHS = L(LRHS)
        KB = L(LB)

        DO I=0,NVT-1
          DWORK (KRHS+I) = DWORK (KRHS+I)*DWORK(KB+I)
        END DO
      
C Release the mass matrix structure again, that's all

        CALL ZDISP (0,LCOLD,'KCOLD ')
        
      ELSE
      
C Ok, we are in the case that we want to generate a real mass
C matrix. Initialize the handles to zero:

        LCOL = 0
        LLD  = 0
        LM   = 0
        NA = 0

C Then create the structure and the matrix:

        ISYMM = 0
        ICLR = 1

        CALL XAP7X(LCOL,LLD,NA,NEQM,TRIA,E011,ISYMM)

        IC = ABS(IPARAM(OICURGS))
        CALL XAB7X(LM,LCOL,LLD,NA,NEQM,
     *             NBLOC,ICLR,TRIA,E011,.FALSE.,COEMGS,
     *             BCONM,KABM,KABMN,IC,ISYMM,ARRDM,0,0D0)
     
C Generate RHS by multiplication of Finite-Difference like LB:

        CALL LAX17 (DWORK(L(LM)),
     *              KWORK(L(LCOL)),KWORK(L(LLD)),NEQM,
     *              DWORK(L(LB)),DWORK(L(LRHS)),1D0,0D0)

C Release the mass matrix again, finish
      
        CALL ZDISP (0,LM,'DM    ')
        CALL ZDISP (0,LCOL,'KCOL  ')
        CALL ZDISP (0,LLD,'KLD   ')
        
      END IF

C Restore the triangulation

C      CALL TRIA2C (TRIOLD)
        
      END

************************************************************************
* Prepare exact right hand side from piecewise linear functions
*
* This routine builds an FE-vector of teh the exact right hand side
* of the static grid deformation. It's an alternative to PRRLG2,
* which only builds the RHS approximatively by multiplying
* with the mass matrix.
*
* The i'th entry of the RHS in the static grid deformation in 
* weak formulation is defined by
*
*   RHS_i = ( 1/f - 1/g , phi_i )
*
* with f being a monitor function, g being the cell distribution and
* phi_i being the test function. This routine now performs a real
* integration to build these entries.
*
* DMON points to the first function f. DMNOLD points to the cell 
* distribution  function g. Both functions consist of the values of 
* the corresponding function in the NVT nodes of the triangulation.
* DRHS is the starting address of the destination vector in DWORK that
* receives the computed right hand side.
*
* The cubature formula is configured by the entry in the IPARAM block.
*
* In:
*   IPARAM : array [1..SZGLSI] of integer = TIGridAdaptionSolver
*            Integer parameter array to configure the solver
*   DPARAM : array [1..SZGLSD] of double = TDGridAdaptionSolver
*            Double precision parameter array to configure the solver
*   TRIA   : array [1..SZTRIA] of integer
*            Triangulation information of the current grid
*   DMON   : array [1..NEQ] of double
*            Monitor function.
*   DMNOLF : array [1..NEQ] of double
*            Old cell distribution
*   DRHS   : array [1..NEQ} of double
*            Right hand side
*   NEQ    : length of vector identified by LRHS
*   ELE    : Element to use for building the right hand side.
*            DMON and DMNOLD are interpreted as coefficient
*            vectors of a FE-function, discretized by ELE.
* Out:
*   The vector specified by DRHS is filled with data.
************************************************************************

      SUBROUTINE PRRLX2 (IPARAM,DPARAM,TRIA,DMON,DMNOLD,DRHS,NEQ,ELE)

      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'cerr.inc'
      INCLUDE 'stria.inc'

      INCLUDE 'cbasictria.inc'

      INCLUDE 'cbasicelem.inc'
      INCLUDE 'ccub.inc'
      INCLUDE 'celem.inc'

      INCLUDE 'sgridadaptsolverparams.inc'
      
C parameter
      
      INTEGER NEQ, TRIA(SZTRIA)
      DOUBLE PRECISION DMON(NEQ),DMNOLD(NEQ),DRHS(NEQ)
      INTEGER IPARAM(SZGLSI)
      DOUBLE PRECISION DPARAM(SZGLSD)

C externals

      EXTERNAL ELE
      INTEGER NDFL
      EXTERNAL NDFL
      
C local variables

      INTEGER KCORVG,IELTYP,IDFL,I,IG,JP,IEL1,IVE,NVE
      DOUBLE PRECISION DJF(2,2),OM,DINTF,XI1,XI2,DF,DG
      DOUBLE PRECISION X,Y,DCOORD(2,4)
      INTEGER KDFG(NNBAS),KDFL(NNBAS)
      INTEGER KVERT, KMID

C     Resolve some handles for quicker access:
      
      KVERT = L(TRIA(OLVERT))
      KMID = L(TRIA(OLMID))
      NVE = TRIA(ONVE)
      
C     We need only function values for the integration:

      DO  I = 1,NNDER
        BDER(I)=.FALSE.
      ENDDO
      BDER(1)=.TRUE.

C     Get the element type as well as number of local degrees of
C     freedom:

      IELTYP=-1
      CALL ELE(0D0,0D0,IELTYP)
      IDFL=NDFL(IELTYP)
      
C     initialise cubature formula
        
      ICUBP = ABS(IPARAM(OICURGS))

      CALL CB2Q(ICUBP)
      IF (IER.NE.0) GOTO 99999

C     Prepare conforming element for cubature. Save the current cubature
C     formula into the COMMON block variable ICUBP - either for now or
C     for later use. Perform a dummy call to the element in the conforming
C     case.
C     This is done for saving arithmetic operations in later calls.

      CALL ELE(0D0,0D0,-2)
        
C     Clear the vector we want to build:

      CALL LCL1(DRHS,NEQ)

C     Loop about the elements, calculate the integral into DINTF

      DINTF = 0D0
      
      DO IEL1 = 1,TRIA(ONEL)
      
C       Calculate the local and global DOF's on our current element.
C       We later have to loop about them...

        CALL NDFGLX(TRIA,IEL1,1,IELTYP,KWORK(KVERT),KWORK(KMID), 
     *              KDFG,KDFL)
        IF (IER.LT.0) GOTO 99999
      
C       Store the coordinates of the corners of the current element in 
C       the COMMON block variables DX/DY (necessary for the element) 
C       as well as in the DCOORD array for the coordinate transformation.

        DO IVE = 1, NVE
          JP=KWORK(KVERT + (IEL1-1)*NNVE+(IVE-1))
          KVE(IVE)=JP
          CALL NDE2XY (JP,TRIA,DCOORD (1,IVE),DCOORD (2,IVE))
          DX(IVE)=DCOORD (1,IVE)
          DY(IVE)=DCOORD (2,IVE)
        END DO
        
C       Initialise auxiliary Jacobian factors DJF for transformation

        CALL QINIJF (DCOORD,DJF)

C       Loop over all cubature points

        DO ICUBP = 1, NCUBP
          
C         Cubature point on the reference element

          XI1=DXI(ICUBP,1)
          XI2=DXI(ICUBP,2)
          
C         Calculate Jacobian matrix of transformation and its
C         determinant. This is necessary for the element routine ELE
C         to properly calculate derivatives. The result is directly
C         written into the element COMMON block variables DJAC and DETJ.
 
          CALL QTRDET (DCOORD,DJF,DJAC,DETJ,XI1,XI2)

C         Evaluate the basis functions in the cubature point

          CALL ELE(XI1,XI2,-3)
          IF (IER.LT.0) GOTO 99999
      
C         Add the values of the basis function to calculate the value
C         of the both functions in the current cubature point.
      
          DF = 0D0
          DG = 0D0
          DO I=1,IDFL
            IG=KDFG(I)
            DF = DF + DMON(IG)*DBAS(KDFL(I),1)
            DG = DG + DMNOLD(IG)*DBAS(KDFL(I),1)
          END DO
          
C         Ok, DF is the value of DMON in the cubature point and
C         DG is the value of DMNOLD there. Now it comes to the
C         integration.
          
C         Calculate the weighting factor for the current cubature point
C         with the help of the Jacobian determinant

          OM = DOMEGA(ICUBP)*DETJ
          
C         Now calculate the value of the right hand side in the current
C         cubature point, weighted by the weight of the cubature formula:

          DINTF = OM * (1D0/DF - 1D0/DG)
C          DINTF = OM
          
C         Remember, the i'th entry of the RHS is defined as:
C
C         r_i = int_Omega (f*phi_i)
C             = sum_T sum_p omega_p * fct(p) * phi_i(p)
C
C         with p=cubature points in the element T and fct(p)=1/f+1/g.
C         phi_i can still be found in the DBAS-array.
C         So now perform a second loop over the basis functions, multiply
C         the calculated value of the function with the basis function
C         and add it to the global degree of freedom to build the RHS.
C         DINTF has already been multiplied with the weight of the
C         integral above.
          
          DO I=1,IDFL
            IG=KDFG(I)
            DRHS(IG) = DRHS(IG) + DINTF*DBAS(KDFL(I),1)
          END DO
          
        END DO
      
      END DO
      
99999 END
      
************************************************************************
* Helper functions for multigrid solver in grid adaption
************************************************************************

************************************************************************
* Callback-routine of multigrid:
* Perform NSM smoothing steps to system A*DX=DB
*
* In:
*   DX     - vector to smooth
*   DB     - RHS of current system
*   DD     - auxilary vector
*   NEQ    - length of vectors
*   NSM    - number of smoothing steps
*   IPARAM  - array [1..SZMGRI] of integer
*             Integer parameter structure for multigrid solver
*   DPARAM  - array [1..SZMGRI] of integer
*             Double precision parameter structure for multigrid solver
*   IDATA   - array [1..*] of integer
*             User defined data array - not used.
*   DDATA   - array [1..*] of double
*             User defined data array - not used.
************************************************************************

      SUBROUTINE YSMG2(DX,DB,DD,NEQ,NSMS,IPARAM, DPARAM,IDATA,DDATA)  
      
      IMPLICIT NONE

      INCLUDE 'cmem.inc'
      INCLUDE 'cout.inc'
      
      INCLUDE 'cbasicmg.inc'
      
      INCLUDE 'stria.inc'

      INCLUDE 'ssolvers.inc'
      INCLUDE 'm020.inc'
      INCLUDE 'sgridadaptstaticsolver.inc'

      INTEGER NEQ, NSMS
      DOUBLE PRECISION DX(NEQ),DB(NEQ),DD(NEQ)
      
      EXTERNAL I000,SM0CG2,YJAX17,YNMFLT
      
      INTEGER IPARAM(SZGSSI)
      DOUBLE PRECISION DPARAM(SZGSSD)
      INTEGER IDATA(*)
      DOUBLE PRECISION DDATA(*)
      
C local variables

      INTEGER KTRA1,KTRA2,NIT,ILEV, LA, LCOL, LLD, LWORK, K
      INTEGER STARTI, STARTD
      DOUBLE PRECISION OMEGA

C Cancel immediately if no smoother is used.
C Should not happen as this is a VERY bad choice for a multigrid
C solver ;)

      IF (IPARAM(OSMTAG).EQ.0) RETURN

C If the matrices are resorted, sort the vectors to fit the structure of
C the matrices.

      ILEV = IPARAM(OILEV)

      IF (IPARAM(OISRT).GT.0) THEN
        KTRA1=L(IPARAM(OLTR1+ILEV-1))
        KTRA2=L(IPARAM(OLTR2+ILEV-1))
        CALL VECSRT(DX,DD,KWORK(KTRA1),KWORK(KTRA2),NEQ,1)
        CALL LCP1  (DD,DX,NEQ)
        CALL VECSRT(DB,DD,KWORK(KTRA1),KWORK(KTRA2),NEQ,1)
        CALL LCP1  (DD,DB,NEQ)
      ENDIF
      
      NIT = NSMS
      
C We have to work with the smoother sub-structure
C of the multigrid structure:

      STARTI = OSMSI
      STARTD = OSMSD

C Get the damping parameter from the multigrid structure. It can be found
C in the OMGSM-variable of the smoother sub-structure in the 
C multigrid structure, as prepared in the initialization routines:      
      
      OMEGA = DPARAM(OOMGSM)
      
C Save this value in the OMGPC-parameter of the smoother structure to
C allow the callback-routine of the (BiCGStab-)smoother to access it:

      DPARAM (STARTD+OOMGPC-1) = OMEGA
      
C All the following operations are to prepare the smoother-structure
C at IPARAM(OSMSI) / DPARAM(OSMSD) to meet the current configuration
C how the smoothing should be performed.

C We only allow BiCGStab as a smoother here - or no smoother at all. 
C Other smoothers don't behave very well concerning pure
C Neumann problems :(
C
C The basic smoother structure for the smoother was initialized
C in INMGG2, but the matrix handle has not been initialized yet
C (as he haven't had that). Now we have that handle - for the system
C matrix as well as the ILU-matrix, if ILU is used as a preconditioner
C in a BiCGStab smoother.
C 
C So now it's the time to initialize the smoother (and probably the
C preconditioner of the smoother) with the matrix handles
C they should use in the following call.
C
C At first initialize the matrix handle structure with the handle of
C the matrix on the coarse grid:

      LA   = IPARAM(OLA+ILEV-1)
      LCOL = IPARAM(OLCOL+ILEV-1)
      LLD  = IPARAM(OLLD+ILEV-1) 

      IPARAM(STARTI+OMATI-1)   = LA
      IPARAM(STARTI+OMATI+1-1) = LCOL
      IPARAM(STARTI+OMATI+2-1) = LLD
    
C Then initialize the preconditioner, if necessary.
C The type of preconditioner that should be used is defined in the
C preconditioner tag of the multigrid structure, as initialized
C in the initialization routines:

      IF (IPARAM(STARTI+OPCTAG-1).NE.0) THEN

C As Jacobi, SSOR and ILU0 are standard preconditioners, they expect 
C handles of the matrix at index PRECI/PRECD of "their" solver 
C structure. The handles have been calculated above - we can 
C simply copy them into where the callback routines expect them:
C
C For Jacobi and SSOR as preconditioner in the smoother,
C the preconditioner routine needs the matrix on the 
C current level; these have been found out some lines above.
C
C For ILU0 as a preconditioner, the preconditioner callback routine
C nneds the ILU matrix:

        IF (IPARAM(STARTI+OPCTAG-1).EQ.8) THEN
          LA   = IPARAM(OLILU0+ILEV-1)
          LCOL = IPARAM(OLCOL+ILEV-1)
          LLD  = IPARAM(OLLD+ILEV-1) 
        END IF

        IPARAM(STARTI+OPRECI-1) = LA
        IPARAM(STARTI+OPRECI+1-1) = LCOL
        IPARAM(STARTI+OPRECI+2-1) = LLD
        
      END IF

      CALL LCP1(DB,DD,NEQ)
      CALL YJAX17(DX,DD,NEQ,-1D0,1D0,IPARAM(STARTI), DPARAM(STARTD))
      CALL LL21(DD,NEQ,OMEGA)

C Reserve some memory for the BiCGStab smoother:
C This is dangerous at the first glance: The array must be deleted
C afterwards, but we have to be sure, that the offset positions of
C the mmultigrid solver does not get destroyed!
C However, since we allocate a new array at the end of the heap,
C deleting it afterwards will not modify any value in the L-array!

      CALL ZNEW (5*NEQ,1,LWORK,'BICGST')
      K = L(LWORK)

C Call the extended BiCGStab-smoother with the smoother structure,
C which determines its behaviour. 
C As matrix-vector multiplication we can use the standard 
C multiplication routine for matrix structure 7. 
C
C IPARAM/DPARAM point to the solver structure of the MG-solver,
C IPARAM(OSMSI)/DPARAM(OSMSD) to the smoothing structure that
C define the behaviour of our smoother.
C This is the structure we have to pass to the smoother.

      IF (IPARAM(OMCGTRM).GE.3) 
     *  WRITE (MTERM,*) 'Smoothing-Step with BiCGStab.'
     
      CALL II21X(IPARAM(STARTI), DPARAM(STARTD), 
     *     IDATA,DDATA,NEQ,DX,DB,
     *     DWORK(K),DWORK(K+NEQ),DWORK(K+2*NEQ),
     *     DWORK(K+3*NEQ),DWORK(K+4*NEQ),
     *     YJAX17,SM0CG2,YNMFLT, NIT)

C Release the memory again.

      CALL LCP1(DB,DD,NEQ)
      CALL YJAX17(DX,DD,NEQ,-1D0,1D0,IPARAM(STARTI), DPARAM(STARTD))
      CALL LL21(DD,NEQ,OMEGA)

      CALL ZDISP (0,LWORK,'BICGST')

C Sort the vectors back if necessary

      IF (IPARAM(OISRT).GT.0) THEN
        CALL VECSRT(DX,DD,KWORK(KTRA1),KWORK(KTRA2),NEQ,2)
        CALL LCP1  (DD,DX,NEQ)
        CALL VECSRT(DB,DD,KWORK(KTRA1),KWORK(KTRA2),NEQ,2)
        CALL LCP1  (DD,DB,NEQ)
      ENDIF

99999 END

************************************************************************
* Callback-routine of multigrid:
* Performs a matrix vector multiplication of the form
*                   DAX:= A1*(A*DX) + A2*DAX
* on the current level with the corresponding system matrix.
* Current system dimension: NEQ.
*
* This matrix-vector multiplication is only designed for being
* used inside of the multigrid algorithm, since it performs
* resorting of the vectors before multiplication with the
* matrix!
************************************************************************

      SUBROUTINE YAXG2(DX,DAX,NEQ,A1,A2,IPARAM,DPARAM,IDATA,DDATA)  

      IMPLICIT NONE

      INCLUDE 'cmem.inc'
      INCLUDE 'cout.inc'
      
      INCLUDE 'cbasicmg.inc'

      INCLUDE 'stria.inc'

      INCLUDE 'ssolvers.inc'
      INCLUDE 'm020.inc'
      INCLUDE 'sgridadaptstaticsolver.inc'

      INTEGER NEQ
      DOUBLE PRECISION DX(NEQ),DAX(NEQ),A1,A2

      INTEGER IPARAM(SZGSSI)
      DOUBLE PRECISION DPARAM(SZGSSD)
      INTEGER IDATA(*)
      DOUBLE PRECISION DDATA(*)
      
C local variables

      INTEGER KTRA1,KTRA2,KAUX,ILEV,LA,LCOL,LLD

C If the matrices are resorted, sort the vectors to fit the structure of
C the matrices.

      ILEV = IPARAM(OILEV)

      IF (IPARAM(OISRT).GT.0) THEN
        KTRA1=L(IPARAM(OLTR1+ILEV-1))
        KTRA2=L(IPARAM(OLTR2+ILEV-1))
        KAUX =L(IPARAM(OLAUX2+ILEV-1))
        CALL VECSRT(DX,DWORK(KAUX),KWORK(KTRA1),KWORK(KTRA2),NEQ,1)
        CALL LCP1  (DWORK(KAUX),DX,NEQ)
        CALL VECSRT(DAX,DWORK(KAUX),KWORK(KTRA1),KWORK(KTRA2),NEQ,1)
        CALL LCP1  (DWORK(KAUX),DAX,NEQ)
      ENDIF

C Matrix-vector-multiplication

      LA   = IPARAM(OLA+ILEV-1)
      LCOL = IPARAM(OLCOL+ILEV-1)
      LLD  = IPARAM(OLLD+ILEV-1) 
      CALL LAX17(DWORK(L(LA)),KWORK(L(LCOL)),
     *           KWORK(L(LLD)),NEQ,DX,DAX,A1,A2)

C Sort the vectors back if necessary

      IF (IPARAM(OISRT).GT.0) THEN
        CALL VECSRT(DX,DWORK(KAUX),KWORK(KTRA1),KWORK(KTRA2),NEQ,2)
        CALL LCP1  (DWORK(KAUX),DX,NEQ)
        CALL VECSRT(DAX,DWORK(KAUX),KWORK(KTRA1),KWORK(KTRA2),NEQ,2)
        CALL LCP1  (DWORK(KAUX),DAX,NEQ)
      END IF

      END

************************************************************************
* Callback-routine of multigrid:
* Set Dirichlet boundary components of DX to 0.
*
* Because we have a pure Neumann problem here, this subroutine
* is of no use!
************************************************************************

      SUBROUTINE YDBCG2(DX,NEQ,IPARAM, DPARAM,IDATA,DDATA)  

      IMPLICIT NONE
      
      INCLUDE 'stria.inc'

      INCLUDE 'cbasicmg.inc'
      INCLUDE 'ssolvers.inc'
      INCLUDE 'm020.inc'
      INCLUDE 'sgridadaptstaticsolver.inc'
      
      INTEGER IPARAM(SZGSSI)
      DOUBLE PRECISION DPARAM(SZGSSD)
      INTEGER IDATA(*)
      DOUBLE PRECISION DDATA(*)

      INTEGER NEQ
      DOUBLE PRECISION DX(NEQ)
      
      
      END

************************************************************************
* Callback-routine of multigrid:
* Compute on the current level ILEV the solution of A*DX=DB
* with an iterative solver.
*
* This calls either BiCGStab or UMFPACK4 as a coarse grid solver.
* Both types of coarse grid solvers are already prepared for usage
* in the initialization routine, so we can directly call them with
* the prepared parameter blocks.
************************************************************************

      SUBROUTINE YEXG2(DX,DB,DD,NEQ,RHO,ITE,IPARAM,DPARAM,IDATA,DDATA)

      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'cout.inc'
      
      INCLUDE 'cbasicmg.inc'

      INCLUDE 'stria.inc'

      INCLUDE 'ssolvers.inc'
      INCLUDE 'm020.inc'
      INCLUDE 'sgridadaptstaticsolver.inc'

      INTEGER NEQ,KSYS,ITE,STARTI,STARTD
      DOUBLE PRECISION DX(NEQ),DB(NEQ),DD(NEQ),RHO
      EXTERNAL I000,MG0CG2,YJAX17,YNMFLT
      
      INTEGER IPARAM(SZGSSI)
      DOUBLE PRECISION DPARAM(SZGSSD)
      INTEGER IDATA(*)
      DOUBLE PRECISION DDATA(*)
      
C local variables

      INTEGER KTRA1,KTRA2,ILEV,K,LWORK

      ILEV = IPARAM(OILEV)

C If the matrices are resorted, sort the vectors to fit the structure of
C the matrices.

      IF (IPARAM(OISRT).GT.0) THEN
        KTRA1=L(IPARAM(OLTR1+ILEV-1))
        KTRA2=L(IPARAM(OLTR2+ILEV-1))
        CALL VECSRT(DX,DD,KWORK(KTRA1),KWORK(KTRA2),NEQ,1)
        CALL LCP1  (DD,DX,NEQ)
        CALL VECSRT(DB,DD,KWORK(KTRA1),KWORK(KTRA2),NEQ,1)
        CALL LCP1  (DD,DB,NEQ)
      ENDIF
      
C We have to work with the structure of the coarse grid solver:

      STARTI = OCGSI
      STARTD = OCGSD
      
C Check the tag in the coarse-grid solver structure which coarse grid
C solver should be used:
      
      IF (IPARAM(STARTI+OSLTAG-1).EQ.6) THEN
C BiCGStab
        IF (IPARAM(OMCGTRM).GE.1) 
     *    WRITE (MTERM,'(A)') 'Solving coarse grid with BiCGStab.'
     
C Reserve some memory for the BiCGStab solver.
C This is dangerous at the first glance: The array must be deleted
C afterwards, but we have to be sure, that the offset positions of
C the mmultigrid solver does not get destroyed!
C However, since we allocate a new array at the end of the heap,
C deleting it afterwards will not modify any value in the L-array!

        CALL ZNEW (5*NEQ,1,LWORK,'BICGST')
        K = L(LWORK)

C All solver preparations have already been done in the initialization
C routine. For the BiCGStab-solver a sub-structure was created
C that contains all handles and properties how it should behave.
C So we can directly call the solver to let it work.
C We only have to take care of, that a filtering routine is
C needed for our grid adaption, as we are solving a Neumann-problem
C on the coarsest level. The calling convention for filtering
C routines is the same for the BiCGStab-solver as well as for the
C multigrid solver. So we can use the multigrid filtering routine
C for filtering in BiCGStab as well. The standard Neumann-filter
C ignores the additional parameters that are passed by multigrid,
C so we can use it directly.
C
C For matrix vector multiplication we can use the standard
C structure-7 routine, as the vectors are already correctly
C sorted and the solver structure is initialized with the
C handles for a structure-7 matrix.
C
C IPARAM/DPARAM point to the solver structure of the MG-solver,
C IPARAM(OCGSI)/DPARAM(OCGSD) to that of the coarse grid solver:
     
        CALL II01X(IPARAM(STARTI), DPARAM(STARTD), 
     *             IDATA,DDATA,NEQ,DX,DB,
     *             DWORK(K),DWORK(K+NEQ),DWORK(K+2*NEQ),
     *             DWORK(K+3*NEQ),DWORK(K+4*NEQ),
     *             YJAX17,MG0CG2,YNMFLT)

C Release the memory again.

        CALL ZDISP (0,LWORK,'BICGST')

        ITE = IPARAM(OITE)
        RHO = DPARAM(ORHO)

      ELSE IF(IPARAM(STARTI+OSLTAG-1).EQ.13) THEN
      
        IF (IPARAM(OMCGTRM).GE.1) 
     *    WRITE (MTERM,'(A)') 'Solving coarse grid with UMFPACK4.'

C UMFPACK4 direct solver
C
C The UMFPACK4 solver needs the CNTRU4/INFOU4 structures we precomputed
C earlier as well as the handle KU4NUM of the factorized matrix.
C Both can be found in the extended parameter blocks of IPARAM as
C all information was precomputed in the initialization/preparation
C routines.
C
C To solve the system, we have to solve the transposed system!
C UMFPACK4 supposes the matrix to be in compressed column format. Our
C storage technique 7 is compressed row format!
C So in fact the matrix we just have factorized is the transposed matrix
C to our origional matrix. Therefore we have to "retranspose" it, i.e.
C we solve the transposed system.

        KSYS = 1
        CALL UMF4SOL(KSYS,DX,DB,IPARAM(OKU4NUM),DPARAM(OCNTRU4),
     *               DPARAM(OINFOU4))
        
C The solution has to be modified to fulfill integral mean value=0

        CALL IMVZER (DX,NEQ)

        ITE = 1
        RHO = 0D0

      END IF
      
C Sort the vectors back if necessary

      IF (IPARAM(OISRT).GT.0) THEN
        CALL VECSRT(DX,DD,KWORK(KTRA1),KWORK(KTRA2),NEQ,2)
        CALL LCP1  (DD,DX,NEQ)
        CALL VECSRT(DB,DD,KWORK(KTRA1),KWORK(KTRA2),NEQ,2)
        CALL LCP1  (DD,DB,NEQ)
      END IF

99999 END

************************************************************************
* MG0CGA - Preconditioning callback routine in coarse grid solver.
*
* In:
*  DG     - array [1..NEQ] of double
*           The vector that is to be preconditioned
*  NEQ    - Length of the vector
*  IPARAM - array [1..SZMGRI] of integer
*           Integer parameter structure of the solver
*  DPARAM - array [1..SZMGRI] of integer
*           Double precision parameter structure of the solver
*  IDATA  - array [1..*] of integer
*           User defined integer array
*  DDATA  - array [1..*] of double
*           User defined double array
*
* Chooses the parameters from the COMMON blocks about the
* preconditioning technique on the coarse grid and performs
* the preconditioning C**(-1)*DG with the preconditioning 
* matrix C**(-1).
*
* This routine is called inside of the coarse grid solver. Therefore
* IPARAM/DPARAM point to the solver structure of the coarse
* grid solver! We have no access to any multigrid information,
* except for that which has been copied by the initialization/
* preparation routines to the coarse grid solver structure - but
* this is all we need here!
************************************************************************

      SUBROUTINE MG0CG2(DG,NEQ,IPARAM, DPARAM,IDATA,DDATA)      

      IMPLICIT NONE

      INCLUDE 'cmem.inc'
      INCLUDE 'cout.inc'

      INCLUDE 'cbasicmg.inc'
      
      INCLUDE 'ssolvers.inc'

      INTEGER NEQ
      DOUBLE PRECISION DG(NEQ)
      
      INTEGER IPARAM(SZSLVI)
      DOUBLE PRECISION DPARAM(SZSLVD)

      INTEGER IDATA(*)
      DOUBLE PRECISION DDATA(*)
      
C local variables

      INTEGER KA,KCOL,KLD,i
      DOUBLE PRECISION OMEGA
      
C The OMEGA-parameter can directly be taken from the solver-structure,
C which points to the structure of the coarse grid solver:
      
      OMEGA = DPARAM(OOMGPC)

C Obtain the handles of the current system matrix from the
C integer structure; this information was copied by the initialization
C routine into the MATI/MATD structure blocks which are reserved
C for matrix information.

      IF(IPARAM(OPCTAG).EQ.1) THEN
C Jacobi preconditioning
        IF (IPARAM(OMSGTRM).GE.3) THEN
          WRITE (MTERM,'(A)') 
     *          'Preconditioning step with Jacobi, MG0CGA.'
        END IF  
 
C Get the information about the current matrix from the matrix-
C information block in the coarse-grid solver structure:

        KA   = L(IPARAM(OMATI))
        KCOL = L(IPARAM(OMATI+1))
        KLD  = L(IPARAM(OMATI+2)) 
        CALL IA117(DWORK(KA),KWORK(KLD),DG,NEQ)
     
      ELSE IF(IPARAM(OPCTAG).EQ.4) THEN
C SSOR      
        IF (IPARAM(OMSGTRM).GE.3) THEN
          WRITE (MTERM,'(A)') 
     *          'Preconditioning step with SSOR, MG0CGA.'
        END IF  
 
C Get the information about the current matrix from the matrix-
C information block in the coarse-grid solver structure:

        KA   = L(IPARAM(OMATI))
        KCOL = L(IPARAM(OMATI+1))
        KLD  = L(IPARAM(OMATI+2)) 
        CALL ID117 (DWORK(KA),KWORK(KCOL),KWORK(KLD),DG,NEQ,
     *              OMEGA) 
     
      ELSE IF(IPARAM(OPCTAG).EQ.8) THEN
C ILU      
        IF (IPARAM(OMSGTRM).GE.3) THEN
          WRITE (MTERM,'(A)') 'Preconditioning step with ILU, MG0CGA.'
        END IF

C The handle of the ILU-matrix can be found in the PRECI part
C of the coarse grid solver structure 

        KA   = L(IPARAM(OPRECI))
        KCOL = L(IPARAM(OPRECI+1))
        KLD  = L(IPARAM(OPRECI+2)) 
        CALL IF117(DWORK(KA),KWORK(KCOL),KWORK(KLD),DG,NEQ) 
        
      ELSE
      
        WRITE (MTERM,'(A)') 'Unknown preconditioner in CGR-solver!?!'
        STOP
        
      ENDIF

99999 END
 
************************************************************************
* SM0CG2 - Preconditioning routine in smoother.
*
* In:
*  DG     - array [1..NEQ] of double
*           The vector that is to be preconditioned
*  NEQ    - Length of the vector
*  IPARAM - array [1..SZMGRI] of integer
*           Integer parameter structure of the solver
*  DPARAM - array [1..SZMGRI] of integer
*           Double precision parameter structure of the solver
*  IDATA  - array [1..*] of integer
*           User defined integer array
*  DDATA  - array [1..*] of double
*           User defined double array

* Chooses the parameters from the COMMON blocks about the
* preconditioning technique in the smoother and performs
* the preconditioning C**(-1)*DG with the preconditioning 
* matrix C**(-1).
*
* This routine is called inside of the smoother. Therefore
* IPARAM/DPARAM point to the smoother structure inside of the MG-
* structure, mot to the MG-structure itself.
* We have no access to any multigrid information,
* except for that which has been copied by the initialization/
* preparation routines to the smoother structure - but
* this is all we need here!
************************************************************************

      SUBROUTINE SM0CG2(DG,NEQ,IPARAM, DPARAM, IDATA,DDATA)      

      IMPLICIT NONE

      INCLUDE 'cmem.inc'
      INCLUDE 'cout.inc'

      INCLUDE 'cbasicmg.inc'
      
      INCLUDE 'ssolvers.inc'

      INTEGER NEQ
      DOUBLE PRECISION DG(NEQ)
      
      INTEGER IPARAM(SZSLVI)
      DOUBLE PRECISION DPARAM(SZSLVD)
      INTEGER IDATA(*)
      DOUBLE PRECISION DDATA(*)
      
C local variables

      INTEGER KA,KCOL,KLD
      DOUBLE PRECISION OMEGA

C The OMEGA-parameter can directly be taken from the "solver-structure",
C which points to the structure of the smoother. The OMGSM-parameter
C from the DAT-file took the following way to here:
C
C  DAT-file  --MG-initialization----------> OMGSM in MG-structure
C  OMGSM     --smoother-callback-routine--> OMGPC in smoother-structure
C  OMGPC     -----------> can be accessed here, DPARAM points to
C                         the smoother-str.

      OMEGA = DPARAM(OOMGPC)

      CALL IMVZER (DG,NEQ)

      IF(IPARAM(OPCTAG).EQ.1) THEN
C Jacobi preconditioning
        IF (IPARAM(OMSGTRM).GE.3) THEN
          WRITE (MTERM,'(A)') 
     *      'Preconditioning step with Jacobi, SM0CGA.'
        END IF  

C Get the information about the current matrix from the matrix-
C information block in the smoother structure:

        KA   = L(IPARAM(OMATI))
        KCOL = L(IPARAM(OMATI+1))
        KLD  = L(IPARAM(OMATI+2)) 
        CALL IA117(DWORK(KA),KWORK(KLD),DG,NEQ)
       
      ELSE IF(IPARAM(OPCTAG).EQ.4) THEN
      
C SSOR      
        IF (IPARAM(OMSGTRM).GE.3) THEN
          WRITE (MTERM,'(A)') 'Preconditioning step with SSOR, SM0CGA.'
        END IF  

C Get the information about the current matrix from the matrix-
C information block in the smoother structure:

        KA   = L(IPARAM(OMATI))
        KCOL = L(IPARAM(OMATI+1))
        KLD  = L(IPARAM(OMATI+2)) 
        CALL ID117 (DWORK(KA),KWORK(KCOL),KWORK(KLD),DG,NEQ,OMEGA) 
     
      ELSE IF(IPARAM(OPCTAG).EQ.8) THEN
      
C ILU0     
        IF (IPARAM(OMSGTRM).GE.3) THEN
          WRITE (MTERM,'(A)') 'Preconditioning step with ILU, SM0CGA.'
        END IF

C Obtain the handles of the current system matrix from the
C integer structure. The preparation routine of the smoother copied
C all necessary information about the preconditioner into the
C preconditioner block of the solver structure for the smoother.

        KA   = L(IPARAM(OPRECI))
        KCOL = L(IPARAM(OPRECI+1))
        KLD  = L(IPARAM(OPRECI+2)) 
        CALL IF117(DWORK(KA),
     *             KWORK(KCOL),KWORK(KLD),DG,NEQ) 

      END IF
      
99999 END

************************************************************************
* Prolongation of vector, E011 element
*
* Extended calling convention. Specialized version of static grid
* adaption.
*
* In:
*   DUC     - array [1..NEQ-coarse] of double
*             coarse grid vector, level ILEV-1, to be prolongated
*   IPARAM  - array [1..SZMGRI] of integer
*             Integer parameter structure for multigrid solver
*   DPARAM  - array [1..SZMGRI] of integer
*             Double precision parameter structure for multigrid solver
*   IDATA   - array [1..*] of integer
*             User defined data array - not used.
*   DDATA   - array [1..*] of double
*             User defined data array - not used.
* Out:
*   DUF     - array [1..NEQ-fine] of double
*             fine grid vector, level ILEV
************************************************************************

      SUBROUTINE YPG2Q1(DUC,DUF,IPARAM,DPARAM,IDATA,DDATA)  
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'cbasicmg.inc'
      
      INCLUDE 'stria.inc'

      INCLUDE 'ssolvers.inc'
      INCLUDE 'm020.inc'
      INCLUDE 'sgridadaptstaticsolver.inc'

      INTEGER IPARAM(SZGSSI)
      DOUBLE PRECISION DPARAM(SZGSSD)
      INTEGER IDATA(*)
      DOUBLE PRECISION DDATA(*)
      
      DOUBLE PRECISION DUF(*),DUC(*)
      
C local variables:

      INTEGER ILEV,KV1,KV2,KA1,KA2,NVT1,NEL1,STARTC,STARTF

      ILEV = IPARAM(OILEV)

C Take the triangulation information from the multigrid structure:

      STARTC = OTRIAS+(ILEV-1-1)*SZTRIA-1
      STARTF = OTRIAS+(ILEV-1)*SZTRIA-1

      KV1=L(IPARAM(STARTC+OLVERT))
      KV2=L(IPARAM(STARTF+OLVERT))
      KA1=L(IPARAM(STARTC+OLADJ))
      KA2=L(IPARAM(STARTF+OLADJ))
      NVT1=IPARAM(STARTC+ONVT)
      NEL1=IPARAM(STARTC+ONEL)

      CALL MP011(DUC,DUF,KWORK(KV1),KWORK(KV2),
     *     KWORK(KA1),KWORK(KA2),NVT1,NEL1)

99999 END

************************************************************************
* Restriction of defect vector, E011 element
*
* Extended calling convention. Specialized version of static grid
* adaption.
*
* In:
*   DDF     - array [1..NEQ-fine] of double
*             fine grid defect vector, level ILEV+1, to be restricted
*   IPARAM  - array [1..SZMGRI] of integer
*             Integer parameter structure for multigrid solver
*   DPARAM  - array [1..SZMGRI] of integer
*             Double precision parameter structure for multigrid solver
*   IDATA   - array [1..*] of integer
*             User defined data array - not used.
*   DDATA   - array [1..*] of double
*             User defined data array - not used.
* Out:
*   DDC     - array [1..NEQ-coarse] of double
*             coarse grid defect vector, level ILEV
************************************************************************

      SUBROUTINE YRG2Q1(DDF,DDC,IPARAM,DPARAM,IDATA,DDATA)  
      
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'cbasicmg.inc'
      
      INCLUDE 'stria.inc'

      INCLUDE 'ssolvers.inc'
      INCLUDE 'm020.inc'
      INCLUDE 'sgridadaptstaticsolver.inc'

      INTEGER IPARAM(SZGSSI)
      DOUBLE PRECISION DPARAM(SZGSSD)
      INTEGER IDATA(*)
      DOUBLE PRECISION DDATA(*)

      DOUBLE PRECISION DDF(*),DDC(*)
      
C local variables:

      INTEGER ILEV,KV1,KV2,KA1,KA2,NVT1,NEL2,STARTC,STARTF

      ILEV = IPARAM(OILEV)

C Take the triangulation information from the multigrid structure:

      STARTC = OTRIAS
      STARTC = OTRIAS+(ILEV-1)*SZTRIA-1
      STARTF = OTRIAS+(ILEV+1-1)*SZTRIA-1

      KV2=L(IPARAM(STARTF+OLVERT))
      KA2=L(IPARAM(STARTF+OLADJ))
      NVT1=IPARAM(STARTC+ONVT)
      NEL2=IPARAM(STARTF+ONEL)

      CALL MR011(DDF,DDC,KWORK(KV2),
     *     KWORK(KA2),NVT1,NEL2)

99999 END


************************************************************************
* Step length control, E011 element
*
* Extended calling convention. Specialized version of static grid
* adaption.
*
* In:
*   DX      - array [1..NEQ] of double
*             Previous solution vector on current level
*   DD      - array [1..NEQ] of double
*             fine correction vector on current level
*   DB      - array [1..NEQ] of double
*             Right hand side vector on current level
*   NEQ     - length of fine-grid vectors
*   IPARAM  - array [1..SZMGRI] of integer
*             Integer parameter structure for multigrid solver
*   DPARAM  - array [1..SZMGRI] of integer
*             Double precision parameter structure for multigrid solver
*   IDATA   - array [1..*] of integer
*             User defined data array - not used.
*   DDATA   - array [1..*] of double
*             User defined data array - not used.
* Out:
*   DALPHA  - Relaxation parameter for coarse grid correction,
*             according to some optimization criterion
************************************************************************

      SUBROUTINE YSG2Q1(DX,DD,DB,NEQ,ALPHA,IPARAM,DPARAM,IDATA,DDATA)  

      IMPLICIT NONE
      
      INTEGER NEQ
      DOUBLE PRECISION ALPHA
      DOUBLE PRECISION DX(NEQ),DD(NEQ),DB(NEQ)
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'stria.inc'
      INCLUDE 'ssolvers.inc'
      INCLUDE 'm020.inc'
      INCLUDE 'sgridadaptstaticsolver.inc'

      INTEGER IPARAM(SZGSSI)
      DOUBLE PRECISION DPARAM(SZGSSD)
      INTEGER IDATA(*)
      DOUBLE PRECISION DDATA(*)

C Q1 doesn't use step length control; it's a conforming element,
C therefor 1D0 is (should be) optimal.

      ALPHA = 1D0

      END
