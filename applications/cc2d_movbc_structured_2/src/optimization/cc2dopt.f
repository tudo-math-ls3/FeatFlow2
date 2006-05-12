************************************************************************
      SUBROUTINE CCOPT
************************************************************************
*
*   Explaination: CC2D optimization module
*
*   Works like CC2D, but is used for optimization.
*   Whether or not optimization is used depends on the setting
*   of the variable IOPTTP. IOPTTP=0 switches off any optimization.
*
************************************************************************

      IMPLICIT NONE

      INCLUDE 'cmem.inc'
      INCLUDE 'cout.inc'
      INCLUDE 'cfiles.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'cns.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmg.inc'
      
      INCLUDE 'cinidat.inc'
      INCLUDE 'ciniopt.inc'

      CHARACTER CDATA*60
      CHARACTER COPTDT*60,CGADAT*60,CGEODT*60,CTMP*60

C local variables

      INTEGER MDATA,MFILE,MSHOW,IWORKG,IWMAXG,IWORKI,IWMAXI
      INTEGER IFILEN,ITFILM
      LOGICAL BNXSTP

      DOUBLE PRECISION TTTSUM,TTT0,TTT1,TTTH,TTLIN,TSPRE,TSSOLV,TSPOS
      DOUBLE PRECISION TTOPRE, TTOPT, TTOPOS, TTOPR1, TTOPO1
      DOUBLE PRECISION TTIT0, TTIT1

C=======================================================================
C     Initialization
C=======================================================================

      TTOPRE = 0D0
      TTOPT  = 0D0
      TTOPOS = 0D0

1091  CALL ZTIME(TTTSUM)

      CALL ZTIME(TTT0)

      CALL ZINIT(NNWORK,'feat.msg','data/CC2D.ERR','data/CC2D.PRT',
     *                             'data/CC2D.SYS','data/CC2D.TRC') 
      CALL ZTIME(TTT1)
      TTLC=TTLC+TTT1-TTT0

      MDATA=79
      
      WRITE (CDATA,'(A)') 'data/cc2d.dat'
      IF (IARGC() .GT. 0) THEN
        CALL GETARG(1,CTMP)
        READ (CTMP,'(A)') CDATA
        WRITE (MTERM,'(A,A60)') 'Using initialization file:   ',CDATA
        WRITE (MTERM,9000) 
      END IF
      
C At the very beginning initialise the output structures.

C Use unit number 62 for file output and open the file as prescribed in
C the DAT file.

      CALL RDOUT (MDATA,CDATA,62,MSHOW)
      MFILE = MFILE1
      
C Read parameters of the DAT file

      CALL RDDAT (MDATA,CDATA,MSHOW)

C Read and set up geometry data

      CGEODT = 'geometry.dat'
      CALL RDGEOM (MDATA,CGEODT)
      CALL INIGEO ()

C Read parameters of the optimizer DAT file
C Use a string variable as some compilers have problems using a
C direct string as parameter!

      WRITE (COPTDT,'(A)') 'BPOS.DAT'
      CALL RDOPT (MDATA,COPTDT,1)
      
C Read parameters of the grid smoother DAT file
C Use a string variable as some compilers have problems using a
C direct string as parameter!

      CGADAT = 'gridadapt.dat'
      CALL RDGSMT (MDATA,CGADAT)
  
C Initialise the optimiser before the solver

      CALL ZTIME(TTT0)

      CALL INIOPT ( (IARGC() .GT. 0), CDATA)
      IF (IER.NE.0) GOTO 99998

      WRITE (MTERM,9000) 

      CALL ZTIME(TTT1)
      TTOPRE=TTOPRE+TTT1-TTT0

C initialise all data and read the DAT file

      CALL ZTIME(TTT0)

C Call INIT1; in case of performing not only a single test,
C don't initialise everything as we are goint to do
C some initialisation on our own later.

      CALL INIT1 (MFILE,IWORKG,IWMAXG,IOPTTP.NE.0)
      IF (IER.NE.0) GOTO 99998
      IWORKI=IWORK
      IWMAXI=IWMAX
      
C Open output channels

      CALL INPWFS
      
C Initialise NS time stepping

      CALL INNSTT
      
      CALL ZTIME(TTT1)
      TTLC=TTLC+TTT1-TTT0

C=======================================================================
C *** Call the (non)steady loop to solve
C=======================================================================

C initialise output variables for nonsteady loop:

      ITFILM=IFINIT-1
      IFILEN=0

C Initialise counter variables in /NSCOUN/; these are incremented 
C in the loop

4444  CALL ZTIME(TTIT0)

      NSUBST = 0
      NNONL  = 0
      NMGU   = 0

C Optimisation preprocessing - in case it's necessary

      CALL ZTIME(TTT0)

      CALL OPTPRE
      IF (IER.NE.0) GOTO 99998

      CALL ZTIME(TTT1)
      TTOPRE=TTOPRE+TTT1-TTT0
      TTOPR1=TTT1-TTT0

C=======================================================================
C User defined analysis    
C Only for testing purposes. To be commmented out for real simulations!
C=======================================================================
      
C      CALL UDAN01      

C=======================================================================
     
C Call the calculation routine to calculate the information of interest.
C This will do the optimisation, call of solver, ...

C Normally this routine is only called once, but in case that the
C optimisation process needs a complete restart, BNXSTP is set to TRUE.
C Then  after doing the postprocessing we have to repeat everything
C again!

      CALL OPTCLC (TSPRE,TSSOLV,TSPOS,BNXSTP)
      IF (IER.NE.0) GOTO 99998

C Call the optimisation postprocess - in case it's necessary

      CALL ZTIME(TTT0)

      CALL OPTPOS
      IF (IER.NE.0) GOTO 99998

      CALL ZTIME(TTT1)
      TTOPOS=TTOPOS+TTT1-TTT0
      TTOPO1=TTT1-TTT0
      
      WRITE (MTERM,*) 
      WRITE (MTERM,'(A,F15.5)') 
     *  'Time, last solver preparation: ',TSPRE
      WRITE (MTERM,'(A,F15.5)') 
     *  'Time, last solver iteration:   ',TSSOLV
      WRITE (MTERM,'(A,F15.5)') 
     *  'Time, last postprocessing:     ',TSPOS
      
      CALL ZTIME(TTIT1)

      WRITE (MTERM,'(A,F15.5)') 
     *  'Time, opt. preprocessing:      ',TTOPR1
      WRITE (MTERM,'(A,F15.5)') 
     *  'Time, opt. postprocessing:     ',TTOPO1
      WRITE (MTERM,'(A,F15.5)') 
     *  'Time, full opt. step:          ',TTIT1-TTIT0
      WRITE (MTERM,*) 
     
C If we have to repeat the calculation, start from the scratch...
     
      IF (BNXSTP) GOTO 4444

C=======================================================================
C     Statistics
C=======================================================================

      IF (MSHOW.GE.2) WRITE (MTERM,*)
      IF (MSHOW.GE.0) WRITE (MFILE,*)
      IF (MSHOW.GE.2) WRITE (MTERM,*) 'STATISTICS :'
      IF (MSHOW.GE.0) WRITE (MFILE,*) 'STATISTICS :'

      IF (MSHOW.GE.2) WRITE(MTERM,*) 'NWORK :      ',NWORK
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'NWORK :      ',NWORK
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'IWORKG:      ',IWORKG
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'IWORKG:      ',IWORKG
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'IWMAXG:      ',IWMAXG
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'IWMAXG:      ',IWMAXG
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'IWORKI:      ',IWORKI
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'IWORKI:      ',IWORKI
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'IWMAXI:      ',IWMAXI
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'IWMAXI:      ',IWMAXI
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'IWORK :      ',IWORK
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'IWORK :      ',IWORK
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'IWMAX :      ',IWMAX
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'IWMAX :      ',IWMAX
      IF (MSHOW.GE.2) WRITE(MTERM,9000)
      IF (MSHOW.GE.0) WRITE(MFILE,9000)
C
      CALL ZTIME(TTT1)
      TTTSUM=TTT1-TTTSUM
      TTTH  =TTGRID+TTPOST+TTADF+TTUPW+TTBDR+TTLC+TTMG
      TTLIN =TTADF+TTUPW+TTBDR+TTLC
C
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'total time : ', TTTSUM
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'total time : ', TTTSUM
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'appr. time : ', TTTH
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'appr. time : ', TTTH
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'grid  time : ', TTGRID
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'grid  time : ', TTGRID
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'post  time : ', TTPOST
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'post  time : ', TTPOST
C
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'lin.  time : ', TTLIN
      IF (MSHOW.GE.0) WRITE(MFILE,*) 'lin.  time : ', TTLIN
      IF (MSHOW.GE.2) WRITE(MTERM,*) '-> mavec time : ', TTADF
      IF (MSHOW.GE.0) WRITE(MFILE,*) '-> mavec time : ', TTADF
      IF (MSHOW.GE.2) WRITE(MTERM,*) '-> konv. time : ', TTUPW
      IF (MSHOW.GE.0) WRITE(MFILE,*) '-> konv. time : ', TTUPW
      IF (MSHOW.GE.2) WRITE(MTERM,*) '-> bdry  time : ', TTBDR
      IF (MSHOW.GE.0) WRITE(MFILE,*) '-> bdry  time : ', TTBDR
      IF (MSHOW.GE.2) WRITE(MTERM,*) '-> LC    time : ', TTLC
      IF (MSHOW.GE.0) WRITE(MFILE,*) '-> LC    time : ', TTLC
C
      IF (MSHOW.GE.2) WRITE(MTERM,*) 'mg    time : ', TTMG
      IF (MSHOW.GE.1) WRITE(MFILE,*) 'mg    time : ', TTMG
C
      IF (MSHOW.GE.2) WRITE(MTERM,*) '#substeps  : ', NSUBST
      IF (MSHOW.GE.0) WRITE(MFILE,*) '#substeps  : ', NSUBST
      IF (MSHOW.GE.2) WRITE(MTERM,*) '#nonlinear : ', NNONL
      IF (MSHOW.GE.1) WRITE(MFILE,*) '#nonlinear : ', NNONL
      IF (MSHOW.GE.2) WRITE(MTERM,*) '#mg        : ', NMGU
      IF (MSHOW.GE.1) WRITE(MFILE,*) '#mg        : ', NMGU
      IF (MSHOW.GE.2) WRITE(MTERM,9000)
      IF (MSHOW.GE.1) WRITE(MFILE,9000)

      IF (MSHOW.GE.2) THEN
      WRITE(MTERM,*)
      WRITE(MTERM,*) ' MULTIGRID COMPONENTS [in percent]:'
      WRITE(MTERM,*) ' smoothing     :', 1.D2*TTS/TTMG
      WRITE(MTERM,*) ' solver        :', 1.D2*TTE/TTMG
      WRITE(MTERM,*) ' defect calc.  :', 1.D2*TTD/TTMG
      WRITE(MTERM,*) ' prolongation  :', 1.D2*TTP/TTMG
      WRITE(MTERM,*) ' restriction   :', 1.D2*TTR/TTMG
      WRITE(MTERM,1)
      ENDIF

      IF (MSHOW.GE.0) THEN
      WRITE(MFILE,*)
      WRITE(MFILE,*) ' MULTIGRID COMPONENTS [in percent]:'
      WRITE(MFILE,*) ' smoothing     :', 1.D2*TTS/TTMG
      WRITE(MFILE,*) ' solver        :', 1.D2*TTE/TTMG
      WRITE(MFILE,*) ' defect calc.  :', 1.D2*TTD/TTMG
      WRITE(MFILE,*) ' prolongation  :', 1.D2*TTP/TTMG
      WRITE(MFILE,*) ' restriction   :', 1.D2*TTR/TTMG
      WRITE(MFILE,1)
      ENDIF

      IF (MSHOW.GE.2) THEN
      WRITE (MTERM,'(A,F15.5)') 
     *  'Total time, opt. preprocessing:  ',TTOPRE
      WRITE (MTERM,'(A,F15.5)') 
     *  'Total time, opt. postprocessing: ',TTOPOS
      WRITE (MTERM,'(A,F15.5)') 
     *  'Total time, last optimization:   ',TTOPT
      WRITE (MTERM,1) 
      END IF

      IF (MSHOW.GE.0) THEN
      WRITE (MFILE,'(A,F15.5)') 
     *  'Total time, opt. preprocessing:  ',TTOPRE
      WRITE (MFILE,'(A,F15.5)') 
     *  'Total time, opt. postprocessing: ',TTOPOS
      WRITE (MFILE,'(A,F15.5)') 
     *  'Total time, last optimization:   ',TTOPT
      WRITE (MFILE,1) 
      END IF

      GOTO 99999
C=======================================================================
C     Error case
C=======================================================================
99998 WRITE(MTERM,*) 'IER', IER
      WRITE(MTERM,*) 'IN SUBROUTINE ',SUB

99999 CLOSE(MFILE)

C Clean up the optimiser

      CALL DONOPT ()

C Clean up the geometry

      CALL DONGEO ()

C Clean up all handles that are allocated in INIT1.
C If this produces an error, check the file done.f that handles
C the deallocation of all arrays!

      CALL DONE1 (0)

D     READ *

9000  FORMAT(79('-'))

      END

************************************************************************
************************************************************************
************************************************************************

************************************************************************
* Initialise optimisation structures
*
* Is called after the initialisation in INIT1
*
* In:
*   BXFN - whether to use CFN as basis for a filename for the
*          output data
*   CFN  - string; if BXFN=true, the filename for the output
*          will be build by: "[CFN].optres". Otherwise CFN is
*          ignored and the output file "BPOS.RES" is opened for
*          writing.
*
* Reads the file BPOS.DAT with information about the optimization.
* Opens the file BPOS.RES with unit no. 67 for writing drag/líft-
* results. If BXFN=true, the file "[CFN].RES" is opened instead of
* BPOS.RES.
************************************************************************

      SUBROUTINE INIOPT (BXFN,CFN)
      
      IMPLICIT NONE

      INCLUDE 'ciniopt.inc'
      
C local variables

      CHARACTER CBFN*60,CFN*60
      LOGICAL BXFN
      INTEGER IFMTS, LSTR
      
      INTEGER STNEWC
      EXTERNAL STNEWC

C We are in the first optimization step

      OPTCNT = 1
      
      IXNUM = 0
      IYNUM = 0
      IRNUM = 0
      
C Open the file where to write the output (drag/lift):

      WRITE (CBFN,'(A)') 'BPOS.RES'
      IF (BXFN) THEN
        LSTR = STNEWC (.TRUE.,CFN)
        CALL STCATC (LSTR,.FALSE.,'.optres')
        CALL STPUT (LSTR, CBFN)
        CALL STDIS (LSTR)
      END IF
      
      IFMTS = 1
      CALL  OF0 (67,CBFN,IFMTS)
      
C Write a headline if we perform more that one iteration

      WRITE (67,'(A)') '#     X-Pos                   '//
     *                   '  Y-Pos                  '//
     *                   'X-Radius                '//
     *                   'Y-Radius                '//
     *                   'Rotation                '//
     *                   'Drag                    '//
     *                   'Lift                    '//
     *                   'abs.Lift                '//
     *                   'Volume                  '//
     *                   'Volume-efficiency       '//
     *                   'Step'
      
C The file will be closed in DONOPT.

      END   

************************************************************************
* Clean up optimisation structures
*
* Is called at the end of the optimisation process to clean up used
* structures.
************************************************************************

      SUBROUTINE DONOPT 
      
      IMPLICIT NONE

      INCLUDE 'cmem.inc'
      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'cns.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmg.inc'
      
      INCLUDE 'cinidat.inc'
      INCLUDE 'ciniopt.inc'
      
C Close the file where we wrote our results to.
      
      CLOSE (67)
      
      END
      
************************************************************************
* Optimisation preprocess
*
* Is called in every loop before the solver is initiated.
************************************************************************
      
      SUBROUTINE OPTPRE
      
      IMPLICIT NONE

      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmgpar.inc'

      INCLUDE 'cinidat.inc'
      INCLUDE 'ciniopt.inc'
      
      END       

************************************************************************
* Optimisation calculation
*
* Calculate the values of interest. For this purpose the Navier-Stokes-
* Solver is called. The information of interest can afterwards be found
* in the COMMON blocks.
*
* MFILE,MSHOW,ITFILM,IFILEN must be provided by the user. The content
* of these variables control the output of the solver to files
* (how much is displayed on screen, number of current output file,...).
* These variables might be changed by the solver in the ongoing
* process.
* In TPRE, TSOLV and TPOS this subroutine returns the time that was
* necessary for preprocessing, solving and postprocessing.
*
* Out:
*  BRCLC - Is set to TRUE if the whole computation must be repeated.
************************************************************************

      SUBROUTINE OPTCLC (TPRE,TSOLV,TPOS,BRCLC)
      
      IMPLICIT NONE

      DOUBLE PRECISION PI
      PARAMETER (PI=3.1415926535897931D0)
      
      INCLUDE 'cmem.inc'
      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'cns.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmg.inc'
      INCLUDE 'cmgpar.inc'
      
      INCLUDE 'cmgtria.inc'

      INCLUDE 'cbasictria.inc'
      INCLUDE 'ctria.inc'

      INCLUDE 'cinidat.inc'
      INCLUDE 'ciniopt.inc'
      INCLUDE 'cinigeometry.inc'
      
      INCLUDE 'sgeometries.inc'
      INCLUDE 'scompgeometries.inc'
      
      INCLUDE 'stria.inc'
      INCLUDE 'dstrings.inc'

C parameters

      DOUBLE PRECISION TPRE, TSOLV, TPOS
      LOGICAL BRCLC

C externals

      DOUBLE PRECISION FBDVOL
      INTEGER GCRSIM,GCRSWV
      EXTERNAL FBDVOL,GCRSIM,GCRSWV
      
      EXTERNAL PRBNSS,PRNSNS
      
C local variables

      DOUBLE PRECISION DINFO(3)
      DOUBLE PRECISION DRVOL,DXTMP,DYTMP,DRTMP,DXFIX,DYFIX,DRFIX,H
      INTEGER I,IINFO(1),J,J2,K,N,U1,U2,U3,S,IDPS
      INTEGER TRIA1(SZTRIA), TRIA2(SZTRIA)
      CHARACTER CFNAME*(40)
      DOUBLE PRECISION FINTOL(10),EPS,HSTART(10)
      INTEGER STEPS
      
C Parameter blocks to pass to the callback routine:
      
      INTEGER IPARAM(64)
      DOUBLE PRECISION DPARAM(64)
      
C Parameter blocks for initial and final coordinate

      DOUBLE PRECISION DSTART(64),DEND(64)
      
C     IPARAM / DPARAM has the following structure:

C     IPARAM(1) = Handle LCORVG to undeformed starting grid
C     IPARAM(2) = Bitfield that configures the design parameter space
C                 Bit0=X-coordinate is changed
C                 Bit1=Y-coordinate is changed
C                 Bit2=rotation is changed
C     IPARAM(3) = What to optimize for
C                 1=minimum drag
C                 2=minimum lift
C
C     DPARAM(1) = TPRE
C     DPARAM(2) = TSOLV
C     DPARAM(3) = TPOS
C     DPARAM(4) = minimum X-coordinate
C     DPARAM(5) = maximum X-coordinate
C     DPARAM(4) = minimum Y-coordinate
C     DPARAM(5) = maximum Y-coordinate
C     DPARAM(4) = minimum rotation
C     DPARAM(5) = maximum rotation

C Let's go...
      
      BRCLC = .FALSE.
      
      TPRE  = 0D0
      TSOLV = 0D0
      TPOS  = 0D0

      IF (IOPTTP.EQ.0) THEN
      
C No optimisation, simply calculate a solution to the
C current configuration.
C We don't need any initialisation (as this was previously done
C by INIT1), but we need postprocessing.
C Perform standard initialisation+postprocessing

        CALL NSCALC (1,18,TPRE,TSOLV,TPOS,DINFO,IINFO)
      
C get the reference volume for computing the volume efficiency

        DRVOL = FBDVOL (0)

C Write out the results for drag/lift to a file

        WRITE (67,'(F8.3,2X,F8.3,3E24.13,4E24.13)') 
     *        DCXPOS,DCYPOS,DCRADX+DCRRLX,DCRADY+DCRRLY,DCROT,
     *        REAL(DINFO(1)),
     *        REAL(DINFO(2)),REAL(DINFO(3)),REAL(DINFO(3)/DRVOL)
      
      ELSE IF (IOPTTP.EQ.1) THEN
      
        PRINT *,'IOPTTP=1 currently not supported.'
        STOP

      ELSE IF (IOPTTP.EQ.2) THEN

C Brute-force test. All X-coordinates, Y-coordinates, rotation angles.
      
C Backup the initial grid

        CALL ZNEW(2*KNVT(NLMAX),1,J,'LCORV2')
        CALL LCP1(DWORK(L(KLCVG(NLMAX))),DWORK(L(J)),2*KNVT(NLMAX))
      
C At first get the reference volume for computing the volume efficiency

        DRVOL = FBDVOL (0)

C IXNUM, IYNUM and IRNUM count the current "position".
C "I" receives the flags how the solver should perform. This
C is dependent of the current position of the obstacle, as grid
C deformation is not always necessary and mostly we can use a
C previous solution as start vector!
        
        IXNUM = 0

C For the first run, initialize I for the call of the solver to        
C -Write out unique GMV-files
C -perform complete grid adaption
C The rest is already initialized

        I = 100

C (We slightly enlarge the upper bound because some Fortran compilers
C have problems with comparing doubles in the following DO-loops)
        
        DO DCXPOS = OBXMIN,OBXMAX+1D-10,OBXSTP
          
          IYNUM = 0
          
          DO DCYPOS = OBYMIN,OBYMAX+1D-10,OBYSTP
            
            IRNUM = 0
            
            DO DCROT = OBRMIN,OBRMAX+1D-10,OBRSTP
            
C Update sin/cos(rot)
          
              DCRSIN = SIN(DCROT*PI/180D0)
              DCRCOS = COS(DCROT*PI/180D0)
              
C Update geometries
            
              CALL UPDGEO()
              
C Restore the original grid

              CALL LCP1(DWORK(L(J)),DWORK(L(KLCVG(NLMAX))),
     *                  2*KNVT(NLMAX))

C Rebuild parameter values of boundary nodes on all levels.
C Necessary because the grid adaption routines need this information.
C Only the coordinates of the vertices are shared between levels,
C the information about the parameter values exist separately!

              DO K=NLMIN,NLMAX
                CALL RXBNPR (DWORK(L(KLCVG(K))),KWORK(L(KLVBD(K))),
     *                   KWORK(L(KLBCT(K))),DWORK(L(KLVBDP(K))),
     *                   NBCT)
              END DO

C Start the solver to compute the values
              
              CALL NSCALC (1,I,TPRE,TSOLV,TPOS,DINFO,IINFO)
              IF (IER.NE.0) RETURN
              
C From the second calculation on, let the initialization routines of
C the solver
C -start with 0 as start vectors
C -rebuild all system specific information
C -perform complete grid adaption
C -write out unique GMV files

              I=103

C Write out the results for drag/lift to a file

              WRITE (67,'(F24.16,2X,F24.16,3E24.13,5E24.13)') 
     *          DCXPOS,DCYPOS,DCRADX+DCRRLX,DCRADY+DCRRLY,DCROT,
     *          REAL(DINFO(1)),REAL(DINFO(2)),REAL(ABS(DINFO(2))),
     *          REAL(DINFO(3)),REAL(DINFO(3)/DRVOL)
              
              OPTCNT = OPTCNT + 1
              
              IRNUM = IRNUM + 1
              
            END DO
            
            IYNUM = IYNUM + 1
            
          END DO
          
C We are at the "end" of a coordinate set

          IXNUM = IXNUM + 1
          
        END DO
      
        CALL ZDISP(0,J,'LCORV2')
      
      ELSE IF (IOPTTP.EQ.3) THEN

C Brute-force test. All X-coordinates with fixed Y-coordinate.
C This is only for visualisation of the grid deformation, as it generates
C a sequence of GMV-files...
      
C Backup the initial grid

        CALL ZNEW(2*KNVT(NLMAX),1,J,'LCORV2')
        CALL LCP1(DWORK(L(KLCVG(NLMAX))),DWORK(L(J)),2*KNVT(NLMAX))
      
C At first get the reference volume for computing the volume efficiency

        DRVOL = FBDVOL (0)

C IXNUM, IYNUM and IRNUM count the current "position".
C "I" receives the flags how the solver should behave. This
C is dependent of the current position of the obstacle, as grid
C deformation is not always necessary and mostly we can use a
C previous solution as start vector!
        
C Start: Write out unique GMV-files, perform complete grid adaption.

        I = 101

        IXNUM = 0
        
C Out X-coordinate here: DCXPOS+SIN(DXTMP).
C We save the basic X-coordinate in DXFIX and use it to compute the 
C "real" DCXPOS!

        DXFIX = DCXPOS
        
C (We slightly enlarge the upper bound because some Fortran compilers
C have problems with comparing doubles in the following DO-loops)
        
        DO DXTMP = OBXMIN,OBXMAX+1D-10,OBXSTP
          
          DCXPOS = DXFIX+SIN(DXTMP*PI/180D0)

C Restore the original grid

          CALL LCP1(DWORK(L(J)),DWORK(L(KLCVG(NLMAX))),2*KNVT(NLMAX))

C Rebuild parameter values of boundary nodes on all levels.
C Necessary because the grid adaption routines need this information.
C Only the coordinates of the vertices are shared between levels,
C the information about the parameter values exist separately!

          DO K=NLMIN,NLMAX
            CALL RXBNPR (DWORK(L(KLCVG(K))),KWORK(L(KLVBD(K))),
     *                   KWORK(L(KLBCT(K))),DWORK(L(KLVBDP(K))),
     *                   NBCT)
          END DO

C Start the solver to compute the values

          CALL NSCALC (1,I,TPRE,TSOLV,TPOS,DINFO,IINFO)
          IF (IER.NE.0) RETURN
              
C Write out the results for drag/lift to a file

          WRITE (67,'(F8.3,2X,F8.3,3E24.13,5E24.13)') 
     *          DCXPOS,DCYPOS,DCRADX+DCRRLX,DCRADY+DCRRLY,DCROT,
     *          REAL(DINFO(1)),REAL(DINFO(2)),REAL(ABS(DINFO(2))),
     *          REAL(DINFO(3)),REAL(DINFO(3)/DRVOL)

          OPTCNT = OPTCNT + 1
              
          IXNUM = IXNUM + 1
          
        END DO
      
        CALL ZDISP(0,J,'LCORV2')
      
      ELSE IF (IOPTTP.EQ.4) THEN

C Brute-force test, flying enterprise in a unit square. 
C
C Only the rotation parameters are treated. X-position, Y-position and
C rotation of the object is depending on the rotation parameters.
C The coordinates/rotation are handled such that the object moves in
C a circle in the unit square.
C
C This is basically for the generation of an animation sequence, i.e.
C for visualisation of the grid deformation, as it generates
C a sequence of GMV-files...
      
C Backup the initial grid

        CALL ZNEW(2*KNVT(NLMAX),1,J,'LCORV2')
        CALL LCP1(DWORK(L(KLCVG(NLMAX))),DWORK(L(J)),2*KNVT(NLMAX))
      
C IRNUM count the current "position".
C "I" receives the flags how the solver should behave. This
C is dependent of the current rotation angle of the object.
        
C Start: Write out unique GMV-files, perform complete grid adaption.

        I = 101

        IXNUM = 0
        
C Final coordinates here: DCXPOS+COS(DXTMP), DCYPOS+SIN(DXTMP)
C Rotation = DXTMP + prescribed rotation in the DAT file.
C We save the basic coordinates/rotation in DxFIX and use it to compute 
C the "real" coordinates/rotation!

        DXFIX = DCXPOS
        DYFIX = DCYPOS
        DRFIX = DCROT
        
C (We slightly enlarge the upper bound because some Fortran compilers
C have problems with comparing doubles in the following DO-loops)
        
        DO DRTMP = OBRMIN,OBRMAX+1D-10,OBRSTP
          
C Initialise position and rotation.
          
          DCXPOS = DXFIX+0.25D0*COS(DRTMP*PI/180D0)
          DCYPOS = DYFIX+0.25D0*SIN(DRTMP*PI/180D0)
          DCROT = DRTMP+DRFIX
          DCRSIN = SIN(DCROT*PI/180D0)
          DCRCOS = COS(DCROT*PI/180D0)
          
C transfer the data to the geometry object if there is one

          CALL UPDGEO ()

C Restore the original grid

          CALL LCP1(DWORK(L(J)),DWORK(L(KLCVG(NLMAX))),2*KNVT(NLMAX))

C Rebuild parameter values of boundary nodes on all levels.
C Necessary because the grid adaption routines need this information.
C Only the coordinates of the vertices are shared between levels,
C the information about the parameter values exist separately!

          DO K=NLMIN,NLMAX
            CALL RXBNPR (DWORK(L(KLCVG(K))),KWORK(L(KLVBD(K))),
     *                   KWORK(L(KLBCT(K))),DWORK(L(KLVBDP(K))),
     *                   NBCT)
          END DO

C Start the solver to compute the values

          CALL NSCALC (1,I,TPRE,TSOLV,TPOS,DINFO,IINFO)
          IF (IER.NE.0) RETURN
              
          OPTCNT = OPTCNT + 1
              
          IXNUM = IXNUM + 1
          
        END DO
      
        CALL ZDISP(0,J,'LCORV2')
        
      ELSE IF (IOPTTP.EQ.5) THEN

C Brute-force test, Rotating screws
C
C Only the rotation parameters are treated. X-position, Y-position and
C rotation of the object is depending on the rotation parameters.
C The rotation specifies the movement of the screws.
C
C This is basically for the generation of an animation sequence, i.e.
C for visualisation of the grid deformation, as it generates
C a sequence of GMV-files...
      
C Backup the initial grid

        CALL ZNEW(2*KNVT(NLMAX),1,J,'LCORV2')
        CALL LCP1(DWORK(L(KLCVG(NLMAX))),DWORK(L(J)),2*KNVT(NLMAX))
      
C IRNUM count the current "position".
C "I" receives the flags how the solver should behave. This
C is dependent of the current rotation angle of the object.
        
C Start: Write out unique GMV-files, perform complete grid adaption.

        I = 101

        IXNUM = 0

C Coordinates do not change, object is fixed.         

C transfer the data to the geometry object if there is one

        CALL UPDGEO ()

C Loop for the "rotation" parameters.
        
C (We slightly enlarge the upper bound because some Fortran compilers
C have problems with comparing doubles in the following DO-loops)
        
        DO DRTMP = OBRMIN,OBRMAX+1D-10,OBRSTP
          
C Initialise the screw-objects with the help of rotation

          DYTMP = SIN(DRTMP*PI/(2*180D0))

C          N = GCRSWV (CGEOM(ICGEOM(OGCIDX)), 0D0,DCRADY, 0D0,  DCRADX,
C     *        DCRADY*0.6, DYTMP*DCRADY/2D0, -DRTMP, 4D0, 
C     *        -DCRADY*0.6, DYTMP*DCRADY/2D0, -DRTMP-180D0, 4D0)
C          N = GCRSWV (CGEOM(ICGEOM(OGCIDX+1)), 0D0,-DCRADY, 0D0, DCRADX,
C     *        DCRADY*0.6, DYTMP*DCRADY/2D0, -DRTMP-180D0, 4D0, 
C     *        -DCRADY*0.6, DYTMP*DCRADY/2D0, -DRTMP, 4D0)

C Direct modification of the screw-parameters - therefore the indices
C used here are a bit complex :)
C The screws start on the heap at DWORK(L(IGEOMS(2,1))) after a header.
C The indices in that array can be found in the corresponding integer-
C array, which is represented by the handle IGEOMS(1,1). That array
C saves the index positions at KWORK(L(IGEOMS(1,1)) + 0..#components-1 )

C Index of the first screw inside of the 2-component object:

          I = KWORK(L(IGEOMS(1,1))+OGCIDX-1)

          N = GCRSWV (DWORK(L(IGEOMS(2,1))+I-1), 0D0,DCRADY, 0D0,  
     *        DCRADX, 0, DCRADY*0.6, DYTMP*DCRADY/2D0, -DRTMP, 4D0, 
     *        -DCRADY*0.6, DYTMP*DCRADY/2D0, -DRTMP-180D0, 4D0)

C and the 2nd screw          
          
          I = KWORK(L(IGEOMS(1,1))+OGCIDX-1+1)

          N = GCRSWV (DWORK(L(IGEOMS(2,1))+I-1), 0D0,-DCRADY, 0D0, 
     *        DCRADX, 0, DCRADY*0.6, DYTMP*DCRADY/2D0, -DRTMP-180D0, 
     *        4D0,-DCRADY*0.6, DYTMP*DCRADY/2D0, -DRTMP, 4D0)
          
C Restore the original grid

          CALL LCP1(DWORK(L(J)),DWORK(L(KLCVG(NLMAX))),2*KNVT(NLMAX))

C Rebuild parameter values of boundary nodes on all levels.
C Necessary because the grid adaption routines need this information.
C Only the coordinates of the vertices are shared between levels,
C the information about the parameter values exist separately!

          DO K=NLMIN,NLMAX
            CALL RXBNPR (DWORK(L(KLCVG(K))),KWORK(L(KLVBD(K))),
     *                   KWORK(L(KLBCT(K))),DWORK(L(KLVBDP(K))),
     *                   NBCT)
          END DO

C Start the solver to compute the values

          CALL NSCALC (1,I,TPRE,TSOLV,TPOS,DINFO,IINFO)
          IF (IER.NE.0) RETURN
              
          OPTCNT = OPTCNT + 1
              
          IXNUM = IXNUM + 1
          
        END DO
      
        CALL ZDISP(0,J,'LCORV2')

      ELSE IF (IOPTTP.EQ.5) THEN

C Brute-force test, Rotating screws
C
C Only the rotation parameters are treated. X-position, Y-position and
C rotation of the object is depending on the rotation parameters.
C The rotation specifies the movement of the screws.
C
C This is basically for the generation of an animation sequence, i.e.
C for visualisation of the grid deformation, as it generates
C a sequence of GMV-files...
      
C Backup the initial grid

        CALL ZNEW(2*KNVT(NLMAX),1,J,'LCORV2')
        CALL LCP1(DWORK(L(KLCVG(NLMAX))),DWORK(L(J)),2*KNVT(NLMAX))
      
C IRNUM count the current "position".
C "I" receives the flags how the solver should behave. This
C is dependent of the current rotation angle of the object.
        
C Start: Write out unique GMV-files, perform complete grid adaption.

        I = 101

        IXNUM = 0

C Coordinates do not change, object is fixed.         

C transfer the data to the geometry object if there is one

        CALL UPDGEO ()

C Loop for the "rotation" parameters.
        
C (We slightly enlarge the upper bound because some Fortran compilers
C have problems with comparing doubles in the following DO-loops)
        
        DO DRTMP = OBRMIN,OBRMAX+1D-10,OBRSTP
          
C Initialise the screw-objects with the help of rotation

          DYTMP = SIN(DRTMP*PI/(2*180D0))

C          N = GCRSWV (CGEOM(ICGEOM(OGCIDX)), 0D0,DCRADY, 0D0,  DCRADX,
C     *        DCRADY*0.6, DYTMP*DCRADY/2D0, -DRTMP, 4D0, 
C     *        -DCRADY*0.6, DYTMP*DCRADY/2D0, -DRTMP-180D0, 4D0)
C          N = GCRSWV (CGEOM(ICGEOM(OGCIDX+1)), 0D0,-DCRADY, 0D0, DCRADX,
C     *        DCRADY*0.6, DYTMP*DCRADY/2D0, -DRTMP-180D0, 4D0, 
C     *        -DCRADY*0.6, DYTMP*DCRADY/2D0, -DRTMP, 4D0)

C Direct modification of the screw-parameters - therefore the indices
C used here are a bit complex :)
C The screws start on the heap at DWORK(L(IGEOMS(2,1))) after a header.
C The indices in that array can be found in the corresponding integer-
C array, which is represented by the handle IGEOMS(1,1). That array
C saves the index positions at KWORK(L(IGEOMS(1,1)) + 0..#components-1 )

C Index of the first screw inside of the 2-component object:

          I = KWORK(L(IGEOMS(1,1))+OGCIDX-1)

          N = GCRSWV (DWORK(L(IGEOMS(2,1))+I-1), 0D0,DCRADY, 0D0,  
     *        DCRADX, 0, DCRADY*0.6, DYTMP*DCRADY/2D0, -DRTMP, 4D0, 
     *        -DCRADY*0.6, DYTMP*DCRADY/2D0, -DRTMP-180D0, 4D0)

C and the 2nd screw          
          
          I = KWORK(L(IGEOMS(1,1))+OGCIDX-1+1)

          N = GCRSWV (DWORK(L(IGEOMS(2,1))+I-1), 0D0,-DCRADY, 0D0, 
     *        DCRADX, 0, DCRADY*0.6, DYTMP*DCRADY/2D0, -DRTMP-180D0, 
     *        4D0,-DCRADY*0.6, DYTMP*DCRADY/2D0, -DRTMP, 4D0)
          
C Restore the original grid

          CALL LCP1(DWORK(L(J)),DWORK(L(KLCVG(NLMAX))),2*KNVT(NLMAX))

C Rebuild parameter values of boundary nodes on all levels.
C Necessary because the grid adaption routines need this information.
C Only the coordinates of the vertices are shared between levels,
C the information about the parameter values exist separately!

          DO K=NLMIN,NLMAX
            CALL RXBNPR (DWORK(L(KLCVG(K))),KWORK(L(KLVBD(K))),
     *                   KWORK(L(KLBCT(K))),DWORK(L(KLVBDP(K))),
     *                   NBCT)
          END DO

C Start the solver to compute the values

          CALL NSCALC (1,I,TPRE,TSOLV,TPOS,DINFO,IINFO)
          IF (IER.NE.0) RETURN
              
          OPTCNT = OPTCNT + 1
              
          IXNUM = IXNUM + 1
          
        END DO
      
        CALL ZDISP(0,J,'LCORV2')
        
      ELSE IF (IOPTTP.EQ.6) THEN

C Brute-force test. All X-coordinates, Y-coordinates, rotation angles,
C mesh deformation, produce sequence of GMV-files with results,
C produce sequence of GMV-files with nodal velocities for analyzing.

C Backup the initial grid

        CALL ZNEW(2*KNVT(NLMAX),1,J,'LCORV2')
        CALL LCP1(DWORK(L(KLCVG(NLMAX))),DWORK(L(J)),2*KNVT(NLMAX))
        
C in the first step we don't have a "last grid"

        J2 = 0
      
C At first get the reference volume for computing the volume efficiency

        DRVOL = FBDVOL (0)

C IXNUM, IYNUM and IRNUM count the current "position".
C "I" receives the flags how the solver should perform. This
C is dependent of the current position of the obstacle, as grid
C deformation is not always necessary and mostly we can use a
C previous solution as start vector!
        
C Start: Write out unique GMV-files, perform complete grid adaption.

        I = 100

        IXNUM = 0
        
C (We slightly enlarge the upper bound because some Fortran compilers
C have problems with comparing doubles in the following DO-loops)
        
        DO DCXPOS = OBXMIN,OBXMAX+1D-10,OBXSTP
          
          IYNUM = 0
          
          DO DCYPOS = OBYMIN,OBYMAX+1D-10,OBYSTP
            
            IRNUM = 0
            
            DO DCROT = OBRMIN,OBRMAX+1D-10,OBRSTP
            
C Update sin/cos(rot)
          
              DCRSIN = SIN(DCROT*PI/180D0)
              DCRCOS = COS(DCROT*PI/180D0)
              
C Update geometries
            
              CALL UPDGEO()
              
C Restore the original grid

              CALL LCP1(DWORK(L(J)),DWORK(L(KLCVG(NLMAX))),
     *                  2*KNVT(NLMAX))

C Rebuild parameter values of boundary nodes on all levels.
C Necessary because the grid adaption routines need this information.
C Only the coordinates of the vertices are shared between levels,
C the information about the parameter values exist separately!

              DO K=NLMIN,NLMAX
                CALL RXBNPR (DWORK(L(KLCVG(K))),KWORK(L(KLVBD(K))),
     *                   KWORK(L(KLBCT(K))),DWORK(L(KLVBDP(K))),
     *                   NBCT)
              END DO

C Start the solver to compute the values

              CALL NSCALC (1,I,TPRE,TSOLV,TPOS,DINFO,IINFO)
              IF (IER.NE.0) RETURN
              
C Normally we continue in the solver with:
C Use solution as start vector, write out unique GMV-files,
C perform partial grid-adaption, rebuild all structures

              I=110

C Write out the results for drag/lift to a file

              WRITE (67,'(F8.3,2X,F8.3,3E24.13,5E24.13)') 
     *          DCXPOS,DCYPOS,DCRADX+DCRRLX,DCRADY+DCRRLY,DCROT,
     *          REAL(DINFO(1)),REAL(DINFO(2)),REAL(ABS(DINFO(2))),
     *          REAL(DINFO(3)),REAL(DINFO(3)/DRVOL)
              
C Do we have a last grid? If not, we now have a last grid...

              IF (J2.EQ.0) THEN
C Backup the current (deformed) grid
                CALL ZNEW(2*KNVT(NLMAX),1,J2,'LCORV3')
                CALL LCP1(DWORK(L(KLCVG(NLMAX))),DWORK(L(J2)),
     *                    2*KNVT(NLMAX))
                CALL C2TRIA(TRIA1)
                TRIA1(OLCORVG) = J2
              ELSE
C Use the last and the current grid to calculate the grid velocity
                CALL C2TRIA(TRIA2)
C Reserve memory for the node velocity components
                CALL ZNEW(KNVT(NLMAX),1,U1,'U1    ')
                CALL ZNEW(KNVT(NLMAX),1,U2,'U2    ')
                CALL ZNEW(KNVT(NLMAX),1,U3,'U3    ')
C Calculate node velocities
                CALL UDAN04 (TRIA1,TRIA2,MIN(OBXSTP,OBYSTP),
     *               DWORK(L(U1)),DWORK(L(U2)),DWORK(L(U3)))
C Create a name for the GMV-file
                S = STNEWC (.FALSE.,'#gmv/GSPEED.gmv.')
                CALL STCATI (S,ILEV,4,.TRUE.)
                CALL STCATC (S,.FALSE.,'.')
                CALL STCATI (S,IXNUM,4,.TRUE.)
                CALL STCATC (S,.FALSE.,'.')
                CALL STCATI (S,IYNUM,4,.TRUE.)
                CALL STCATC (S,.FALSE.,'.')
                CALL STCATI (S,IRNUM,4,.TRUE.)
                CALL STCATC (S,.FALSE.,'.')
                CALL STPUT (S, CFNAME)
                CALL STDIS (S)
C Write the node velocities into a GMV-file
                CALL XGMVSC(68,CFNAME,TRIA2(ONEL),TRIA2(ONVT),
     *               KWORK(L(TRIA2(OLVERT))),DWORK(L(TRIA2(OLCORVG))),
     *               DWORK(L(U1)),DWORK(L(U2)),DWORK(L(U3)),
     *               DBLE(OPTCNT))
C Release memory
                CALL ZDISP(0,U3,'U3    ')
                CALL ZDISP(0,U2,'U2    ')
                CALL ZDISP(0,U1,'U1    ')
C Backup the current grid for the next step
                CALL LCP1(DWORK(L(TRIA2(OLCORVG))),
     *                    DWORK(L(TRIA1(OLCORVG))),
     *                    2*TRIA1(ONVT))
              END IF
              
              OPTCNT = OPTCNT + 1
              
              IRNUM = IRNUM + 1
              
            END DO
            
            IYNUM = IYNUM + 1
            
          END DO
          
C We are at the "end" of a coordinate set; the next solver run must
C start with a 0-vector

          I=111
          
          IXNUM = IXNUM + 1
          
        END DO
      
C Release memory

        CALL ZDISP (0,J2,'LCORV2') 
      
      ELSE IF (IOPTTP.EQ.7) THEN

C Brute-force test. All X-coordinates, Y-coordinates, rotation angles,
C mesh deformation, produce sequence of GMV-files with results,
C produce sequence of GMV-files with nodal velocities for analyzing.

C Backup the initial grid

        CALL ZNEW(2*KNVT(NLMAX),1,J,'LCORV2')
        CALL LCP1(DWORK(L(KLCVG(NLMAX))),DWORK(L(J)),2*KNVT(NLMAX))
        
C in the first step we don't have a "last grid"

        J2 = 0
      
C At first get the reference volume for computing the volume efficiency

        DRVOL = FBDVOL (0)

C IXNUM, IYNUM and IRNUM count the current "position".
C "I" receives the flags how the solver should perform. This
C is dependent of the current position of the obstacle, as grid
C deformation is not always necessary and mostly we can use a
C previous solution as start vector!
        
C Start: Write out unique GMV-files, perform complete grid adaption.

        IXNUM = 0
        
        DCXPOS = 0.2
        DCYPOS = 0.2

C Update sin/cos(rot)
          
        DCRSIN = SIN(DCROT*PI/180D0)
        DCRCOS = COS(DCROT*PI/180D0)
              
C Update geometries
            
        CALL UPDGEO()
        
C Start the solver

        OPTCNT = 1
        I = 110       
        CALL NSCALC (1,I,TPRE,TSOLV,TPOS,DINFO,IINFO)
        IF (IER.NE.0) RETURN

C Move the circle

        DCXPOS = 1.5

C Update geometries
            
        CALL UPDGEO()
        
C Start the solver

        OPTCNT = 1
        I = 110       
        CALL NSCALC (1,I,TPRE,TSOLV,TPOS,DINFO,IINFO)
        IF (IER.NE.0) RETURN
        
C Release memory

        CALL ZDISP (0,J2,'LCORV2') 

      ELSE IF (IOPTTP.EQ.8) THEN

C Brute-force test. All X-coordinates, Y-coordinates, rotation angles.
C Don't reset the grid during the computation.
      
C Backup the initial grid

        CALL ZNEW(2*KNVT(NLMAX),1,J,'LCORV2')
        CALL LCP1(DWORK(L(KLCVG(NLMAX))),DWORK(L(J)),2*KNVT(NLMAX))
      
C At first get the reference volume for computing the volume efficiency

        DRVOL = FBDVOL (0)

C IXNUM, IYNUM and IRNUM count the current "position".
C "I" receives the flags how the solver should perform. This
C is dependent of the current position of the obstacle, as grid
C deformation is not always necessary and mostly we can use a
C previous solution as start vector!
        
        IXNUM = 0

C For the first run, initialize I for the call of the solver to        
C -Write out unique GMV-files
C -perform complete grid adaption
C The rest is already initialized

        I = 100

C (We slightly enlarge the upper bound because some Fortran compilers
C have problems with comparing doubles in the following DO-loops)
        
        DO DCXPOS = OBXMIN,OBXMAX+1D-10,OBXSTP
          
          IYNUM = 0
          
          DO DCYPOS = OBYMIN,OBYMAX+1D-10,OBYSTP
            
            IRNUM = 0
            
            DO DCROT = OBRMIN,OBRMAX+1D-10,OBRSTP
            
C Update sin/cos(rot)
          
              DCRSIN = SIN(DCROT*PI/180D0)
              DCRCOS = COS(DCROT*PI/180D0)
              
C Update geometries
            
              CALL UPDGEO()
              
C Don't restore the original grid in this type of calculation!
C              CALL LCP1(DWORK(L(J)),DWORK(L(KLCVG(NLMAX))),
C     *                  2*KNVT(NLMAX))

C Rebuild parameter values of boundary nodes on all levels.
C Necessary because the grid adaption routines need this information.
C Only the coordinates of the vertices are shared between levels,
C the information about the parameter values exist separately!

              DO K=NLMIN,NLMAX
                CALL RXBNPR (DWORK(L(KLCVG(K))),KWORK(L(KLVBD(K))),
     *                   KWORK(L(KLBCT(K))),DWORK(L(KLVBDP(K))),
     *                   NBCT)
              END DO

C Start the solver to compute the values
              
              CALL NSCALC (1,I,TPRE,TSOLV,TPOS,DINFO,IINFO)
              IF (IER.NE.0) RETURN
              
C From the second calculation on, let the initialization routines of
C the solver
C -start with 0 as start vectors
C -rebuild all system specific information
C -perform complete grid adaption
C -write out unique GMV files

              I=103

C Write out the results for drag/lift to a file

              WRITE (67,'(F8.3,2X,F8.3,3E24.13,5E24.13)') 
     *          DCXPOS,DCYPOS,DCRADX+DCRRLX,DCRADY+DCRRLY,DCROT,
     *          REAL(DINFO(1)),REAL(DINFO(2)),REAL(ABS(DINFO(2))),
     *          REAL(DINFO(3)),REAL(DINFO(3)/DRVOL)
              
              OPTCNT = OPTCNT + 1
              
              IRNUM = IRNUM + 1
              
            END DO
            
            IYNUM = IYNUM + 1
            
          END DO
          
C We are at the "end" of a coordinate set

          IXNUM = IXNUM + 1
          
        END DO
      
        CALL ZDISP(0,J,'LCORV2')
      
      ELSE IF (IOPTTP.EQ.100) THEN
      
C       Optimization algorithm 00: Compass search
C
C       First make a backup of the undeformed grid: to use it 
C       as reference:

        CALL ZNEW(2*KNVT(NLMAX),1,J,'LCORV2')
        CALL LCP1(DWORK(L(KLCVG(NLMAX))),DWORK(L(J)),2*KNVT(NLMAX))

C       Build the IPARAM/DPARAM parameter blocks for the 
C       callback routine:

        DPARAM(1) = 0D0
        DPARAM(2) = 0D0
        DPARAM(3) = 0D0
        DPARAM(4) = OBXMIN
        DPARAM(5) = OBXMAX
        DPARAM(6) = OBYMIN
        DPARAM(7) = OBYMAX
        DPARAM(8) = OBRMIN
        DPARAM(9) = OBRMAX

        IPARAM(1) = J
        
C       Search for minimal lift

        IPARAM(3) = 2
        
C       Determine what to optimize for, determine the dimension
C       of the design parameter space and determine the starting
C       position

        IPARAM(2) = 0
        IDPS = 0        
        
        IF (OBXMIN.LT.OBXMAX) THEN
          IPARAM(2) = OR(IPARAM(2),1)
          IDPS = IDPS+1
          HSTART(IDPS) = OBXSTP
          DSTART(IDPS) = OBXINI
        END IF

        IF (OBYMIN.LT.OBYMAX) THEN
          IPARAM(2) = OR(IPARAM(2),2)
          IDPS = IDPS+1
          HSTART(IDPS) = OBYSTP
          DSTART(IDPS) = OBYINI
        END IF
        
        IF (OBRMIN.LT.OBRMAX) THEN
          IPARAM(2) = OR(IPARAM(2),4)
          IDPS = IDPS+1
          HSTART(IDPS) = OBRSTP
          DSTART(IDPS) = OBRINI
        END IF
        
        EPS = DOPEPS

        IF (IPARAM(2).NE.0) THEN
        
C         Call the optimizer:

          CALL OPTA00 (IDPS,DSTART,HSTART,EPS,OITMAX,
     *                 IPARAM,DPARAM,64,PRBNSS,PRNSNS,
     *                 DEND,FINTOL,STEPS)
        
        END IF

C       Get results from DPARAM:

        TPRE  = DPARAM(1) 
        TSOLV = DPARAM(2) 
        TPOS  = DPARAM(3)

C       Restore grid:

        CALL LCP1(DWORK(L(J)),DWORK(L(KLCVG(NLMAX))),
     *            2*KNVT(NLMAX))

C       Rebuild parameter values of boundary nodes on all levels.

        DO K=NLMIN,NLMAX
          CALL RXBNPR (DWORK(L(KLCVG(K))),KWORK(L(KLVBD(K))),
     *                 KWORK(L(KLBCT(K))),DWORK(L(KLVBDP(K))),
     *                 NBCT)
        END DO

C       Release memory, finish

        CALL ZDISP (0,J,'LCORV2') 
      
      ELSE IF (IOPTTP.EQ.101) THEN
      
C       Optimization algorithm 01: Nelder Mead
C
C       First make a backup of the undeformed grid: to use it 
C       as reference:

        CALL ZNEW(2*KNVT(NLMAX),1,J,'LCORV2')
        CALL LCP1(DWORK(L(KLCVG(NLMAX))),DWORK(L(J)),2*KNVT(NLMAX))

C       Build the IPARAM/DPARAM parameter blocks for the 
C       callback routine:

        DPARAM(1) = 0D0
        DPARAM(2) = 0D0
        DPARAM(3) = 0D0
        DPARAM(4) = OBXMIN
        DPARAM(5) = OBXMAX
        DPARAM(6) = OBYMIN
        DPARAM(7) = OBYMAX
        DPARAM(8) = OBRMIN
        DPARAM(9) = OBRMAX

        IPARAM(1) = J
        
C       Search for minimal lift

        IPARAM(3) = 2
        
C       Determine what to optimize for, determine the dimension
C       of the design parameter space and determine the starting
C       position

        IPARAM(2) = 0
        IDPS = 0        
        
        IF (OBXMIN.LT.OBXMAX) THEN
          IPARAM(2) = OR(IPARAM(2),1)
          IDPS = IDPS+1
          HSTART(IDPS) = OBXSTP
          DSTART(IDPS) = OBXINI
        END IF

        IF (OBYMIN.LT.OBYMAX) THEN
          IPARAM(2) = OR(IPARAM(2),2)
          IDPS = IDPS+1
          HSTART(IDPS) = OBYSTP
          DSTART(IDPS) = OBYINI
        END IF
        
        IF (OBRMIN.LT.OBRMAX) THEN
          IPARAM(2) = OR(IPARAM(2),4)
          IDPS = IDPS+1
          HSTART(IDPS) = OBRSTP
          DSTART(IDPS) = OBRINI
        END IF
        
        EPS = DOPEPS

        IF (IPARAM(2).NE.0) THEN
        
C         Call the optimizer:

          CALL OPTA01 (IDPS,DSTART,HSTART,EPS,OITMAX,
     *                 0.5D0, 0.5D0, 2.0D0,
     *                 IPARAM,DPARAM,64,PRBNSS,PRNSNS,
     *                 DEND,FINTOL(1),STEPS)    
        END IF

C       Get results from DPARAM:

        TPRE  = DPARAM(1) 
        TSOLV = DPARAM(2) 
        TPOS  = DPARAM(3)

C       Restore grid:

        CALL LCP1(DWORK(L(J)),DWORK(L(KLCVG(NLMAX))),
     *            2*KNVT(NLMAX))

C       Rebuild parameter values of boundary nodes on all levels.

        DO K=NLMIN,NLMAX
          CALL RXBNPR (DWORK(L(KLCVG(K))),KWORK(L(KLVBD(K))),
     *                 KWORK(L(KLBCT(K))),DWORK(L(KLVBDP(K))),
     *                 NBCT)
        END DO

C       Release memory, finish

        CALL ZDISP (0,J,'LCORV2') 
      
      ELSE IF (IOPTTP.EQ.102) THEN
      
C       Optimization algorithm 02: Line-search with finite-difference-
C                                  like approximate gradient
C
C       First make a backup of the undeformed grid: to use it 
C       as reference:

        CALL ZNEW(2*KNVT(NLMAX),1,J,'LCORV2')
        CALL LCP1(DWORK(L(KLCVG(NLMAX))),DWORK(L(J)),2*KNVT(NLMAX))

C       Build the IPARAM/DPARAM parameter blocks for the 
C       callback routine:

        DPARAM(1) = 0D0
        DPARAM(2) = 0D0
        DPARAM(3) = 0D0
        DPARAM(4) = OBXMIN
        DPARAM(5) = OBXMAX
        DPARAM(6) = OBYMIN
        DPARAM(7) = OBYMAX
        DPARAM(8) = OBRMIN
        DPARAM(9) = OBRMAX

        IPARAM(1) = J
        
C       Search for minimal lift

        IPARAM(3) = 2
        
C       Determine what to optimize for, determine the dimension
C       of the design parameter space and determine the starting
C       position

        IPARAM(2) = 0
        IDPS = 0        
        
        IF (OBXMIN.LT.OBXMAX) THEN
          IPARAM(2) = OR(IPARAM(2),1)
          IDPS = IDPS+1
          HSTART(IDPS) = OBXSTP
          DSTART(IDPS) = OBXINI
        END IF

        IF (OBYMIN.LT.OBYMAX) THEN
          IPARAM(2) = OR(IPARAM(2),2)
          IDPS = IDPS+1
          HSTART(IDPS) = OBYSTP
          DSTART(IDPS) = OBYINI
        END IF
        
        IF (OBRMIN.LT.OBRMAX) THEN
          IPARAM(2) = OR(IPARAM(2),4)
          IDPS = IDPS+1
          HSTART(IDPS) = OBRSTP
          DSTART(IDPS) = OBRINI
        END IF
        
        EPS = DOPEPS

        IF (IPARAM(2).NE.0) THEN
        
C         Call the optimizer:

          CALL OPTA02 (IDPS,DSTART,HSTART,DOPINS,EPS,OITMAX,
     *                 IPARAM,DPARAM,64,PRBNSS,PRNSNS,
     *                 DEND,FINTOL,STEPS)    
        END IF

C       Get results from DPARAM:

        TPRE  = DPARAM(1) 
        TSOLV = DPARAM(2) 
        TPOS  = DPARAM(3)

C       Restore grid:

        CALL LCP1(DWORK(L(J)),DWORK(L(KLCVG(NLMAX))),
     *            2*KNVT(NLMAX))

C       Rebuild parameter values of boundary nodes on all levels.

        DO K=NLMIN,NLMAX
          CALL RXBNPR (DWORK(L(KLCVG(K))),KWORK(L(KLVBD(K))),
     *                 KWORK(L(KLBCT(K))),DWORK(L(KLVBDP(K))),
     *                 NBCT)
        END DO

C       Release memory, finish

        CALL ZDISP (0,J,'LCORV2') 

      END IF

99998 END

************************************************************************
* Probe Navier Stokes solver
*
* This routine is the probing callback routine to calculate a
* value of interest with the Navier Stokes solver.
************************************************************************

      SUBROUTINE PRBNSS (IDIMS,LCOORDS,NRES,LRES,IINFO,DINFO)
      
      IMPLICIT NONE
      
      INCLUDE 'cinigeometry.inc'
      INCLUDE 'cmem.inc'
      INCLUDE 'cout.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmgtria.inc'
      INCLUDE 'cmgpar.inc'
      INCLUDE 'cmgfld.inc'
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'ctria.inc'
      
C parameters

      INTEGER IDIMS, IINFO(*),NRES,LRES
      DOUBLE PRECISION DINFO(*)
      INTEGER LCOORDS

C externals:

      DOUBLE PRECISION FBDVOL
      EXTERNAL FBDVOL

      DOUBLE PRECISION PI
      PARAMETER (PI=3.1415926535897931D0)

C OPTCLC defines the IINFO/DINFO parameter blocks the following way:

C     IPARAM / DPARAM has the following structure:

C     IINFO(1) = Handle LCORVG to undeformed starting grid
C     IINFO(2) = Bitfield that configures the design parameter space
C                 Bit0=X-coordinate is changed
C                 Bit1=Y-coordinate is changed
C                 Bit2=rotation is changed
C     IINFO(3) = What to optimize for
C                 1=minimum drag
C                 2=minimum lift

C     DINFO(1) = TPRE
C     DINFO(2) = TSOLV
C     DINFO(3) = TPOS
C     DINFO(4) = minimum X-coordinate
C     DINFO(5) = maximum X-coordinate
C     DINFO(6) = minimum Y-coordinate
C     DINFO(7) = maximum Y-coordinate
C     DINFO(8) = minimum rotation
C     DINFO(9) = maximum rotation
      
C local variables

      INTEGER I,K
      INTEGER IDATA(64)
      DOUBLE PRECISION RES(64)

C At first we have to define DCXPOS/DCYPOS/DCROT to the new configuration.
C The basic parameter setting can be found in DCOORDS. The Bitfield in 
C IINFO(2) decides which information can be found in DCOORDS.

      I = 0
      
      IF (AND(IINFO(2),1).NE.0) THEN
        I = I+1
        DCXPOS = DWORK(L(LCOORDS)+I-1)
      ELSE
        DCXPOS = DINFO(4)
      END IF

      IF (AND(IINFO(2),2).NE.0) THEN
        I = I+1
        DCYPOS = DWORK(L(LCOORDS)+I-1)
      ELSE
        DCYPOS = DINFO(6)
      END IF

      IF (AND(IINFO(2),4).NE.0) THEN
        I = I+1
        DCROT = DWORK(L(LCOORDS)+I-1)
      ELSE
        DCROT = DINFO(8)
      END IF
              
C Update geometries
            
      DCRSIN = SIN(DCROT*PI/180D0)
      DCRCOS = COS(DCROT*PI/180D0)
      CALL UPDGEO()
              
C If we are outside of the allowed domain, cancel immediately

      DWORK(L(LRES)) = 1D99

      IF ( (DCXPOS.LT.DINFO(4)-1D-14).OR.(DCXPOS.GT.DINFO(5)+1D-14).OR.
     *     (DCYPOS.LT.DINFO(6)-1D-14).OR.(DCYPOS.GT.DINFO(7)+1D-14).OR.
     *     (DCROT .LT.DINFO(8)-1D-14).OR.(DCROT .GT.DINFO(9)+1D-14) ) 
     *     RETURN
              
C Restore the original grid

      CALL LCP1(DWORK(L(IINFO(1))),DWORK(L(KLCVG(NLMAX))),
     *          2*KNVT(NLMAX))

C Rebuild parameter values of boundary nodes on all levels.
C Necessary because the grid adaption routines need this information.
C Only the coordinates of the vertices are shared between levels,
C the information about the parameter values exist separately!

      DO K=NLMIN,NLMAX
        CALL RXBNPR (DWORK(L(KLCVG(K))),KWORK(L(KLVBD(K))),
     *               KWORK(L(KLBCT(K))),DWORK(L(KLVBDP(K))),
     *               NBCT)
      END DO

C Start the solver to compute the values into RES

      CALL NSCALC (1,103,DINFO(1),DINFO(2),DINFO(3),RES,IDATA)
      
C Get the reference volume for computing the volume efficiency,
C put it into RES(5)

      RES(5) = FBDVOL (0)

C Copy results to result array for later use

      CALL LCP1(RES,DWORK(L(LRES)+1),3)
      
C Return result value in first entry of RES:

      IF (IINFO(3).EQ.1) DWORK(L(LRES)) = RES(1)
      IF (IINFO(3).EQ.2) DWORK(L(LRES)) = ABS(RES(2))
      
      WRITE (MTERM,*) 'Current probing value: ',DWORK(L(LRES))
      
      END

************************************************************************
* Probing of Navier Stokes solver: Next iterate
*
* This callback routine is called by optimization subroutines
* if the optimization routine processes to a next step.
************************************************************************

      SUBROUTINE PRNSNS (IDIMS,DCOORDS,DCOOPT,IINFO,DINFO,
     *                   PITER,STEP,IACPT,NRES,DRES,DRSOPT)
      
      IMPLICIT NONE
      
      INCLUDE 'cinigeometry.inc'
      INCLUDE 'cmem.inc'
      INCLUDE 'cout.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmgtria.inc'
      INCLUDE 'cmgpar.inc'
      INCLUDE 'cmgfld.inc'
      
      INCLUDE 'cbasictria.inc'
      INCLUDE 'ctria.inc'
      
      INCLUDE 'ciniopt.inc'
      
      INCLUDE 'dstrings.inc'

C parameters

      INTEGER IDIMS, IINFO(*), STEP, NRES,PITER,IACPT
      DOUBLE PRECISION DCOORDS(IDIMS),DCOOPT(IDIMS),DINFO(*)
      DOUBLE PRECISION DRES(*),DRSOPT(*)

C local variables

      INTEGER I
      DOUBLE PRECISION DCX,DCY,DCR,RES(64),D

      INTEGER ITFILM,IFILEN,IFPOST,MSHOW
      DOUBLE PRECISION DDRAG,DLIFT
      CHARACTER*128 CFNAME
      
      DOUBLE PRECISION UE
      EXTERNAL UE
      
C At first we have to define DCXPOS/DCYPOS/DCROT to the new configuration.
C The basic parameter setting can be found in DCOORDS. The Bitfield in 
C IINFO(2) decides which information can be found in DCOORDS.

      I = 0
      
      IF (AND(IINFO(2),1).NE.0) THEN
        I = I+1
        DCX = DCOORDS(I)
      ELSE
        DCX = DINFO(4)
      END IF

      IF (AND(IINFO(2),2).NE.0) THEN
        I = I+1
        DCY = DCOORDS(I)
      ELSE
        DCY = DINFO(6)
      END IF

      IF (AND(IINFO(2),4).NE.0) THEN
        I = I+1
        DCR = DCOORDS(I)
      ELSE
        DCR = DINFO(8)
      END IF

C Copy result array to RES for easier access

      CALL LCP1(DRES,RES,3)
      
      OPTCNT = STEP

C Write results to file

      CALL LCP1(DRES,RES,5)

      IF (IACPT.NE.0) THEN

        D=RES(5)
        IF (D.EQ.0D0) D = 1D0
        WRITE (67,'(F24.16,2X,F24.16,3E24.13,5E24.13,I24)') 
     *          DCX,DCY,DCRADX,DCRADY,DCR,
     *          REAL(RES(2)),REAL(RES(3)),REAL(ABS(RES(3))),
     *          REAL(RES(4)),REAL(RES(4)/D),PITER
     *          

C Build a filename for the GMV-output:

        I = STNEWC (.FALSE.,'gmv/OPT.gmv.')
        CALL STCATI (I,STEP,4,.TRUE.)
        CALL STPUT (I, CFNAME)

        IF (MT.GE.1) THEN
          WRITE (MTERM,'(A$)') 'Writing GMV-file: '
          CALL STOUT (I,MTERM,.TRUE.)
        END IF

        CALL STDIS (I)

C Perform simple postprocessing

        MSHOW = 2
        IFILEN=0
        IFPOST=0
        ITFILM=0
        I=1

        CALL FPOST(IFPOST,IFILEN,ITFILM,UE,MSHOW,DDRAG,DLIFT,I,CFNAME)

      END IF

C Print out current position

      IF (MT.GE.1) THEN
        WRITE (MTERM,'(A)') 
     *       '=================================================='
        WRITE (MTERM,'(A,I10)') 
     *     'Current iterate:           STEP   = ',STEP
        WRITE (MTERM,'(A,I10)') 
     *     'Current optimization it.:  ITER   = ',PITER
        WRITE (MTERM,'(A,I10)') 
     *     'Type of obstacle:          ICIRTP = ',ICIRTP
        WRITE (MTERM,'(A,F10.4)') 
     *     'Obstacle at position:      DCXPOS = ',DCXPOS
        WRITE (MTERM,'(A,F10.4)') 
     *     '                           DCYPOS = ',DCYPOS
        WRITE (MTERM,'(A,F10.4)') 
     *     'X-Radius of obstacle:      DCRADX = ',DCRADX
        WRITE (MTERM,'(A,F10.4)') 
     *     'Y-Radius of obstacle:      DCRADY = ',DCRADY
        WRITE (MTERM,'(A,F10.4)') 
     *     'Rotation of obstacle:      DCROT  = ',DCROT
        IF (IACPT.EQ.0) THEN
          WRITE (MTERM,'(A,E24.16)') 
     *     'Calculated probing value:  ',DRES(1)
          WRITE (MTERM,'(A,E24.16)') 
     *     'Current best pr. value:    ',DRSOPT(1)
        ELSE
          WRITE (MTERM,'(A,E24.16)') 
     *     'Calculated probing value:  ',DRES(1)
          WRITE (MTERM,'(A,E24.16)') 
     *     'Previous best pr. value:   ',DRSOPT(1)
        END IF
        WRITE (MTERM,'(A,F10.4)') 
     *       '=================================================='
        WRITE (MTERM,'(A)') 
      END IF

      END 

************************************************************************
* Optimisation postprocess
*
* Write out the results of the last calculation.
************************************************************************
      
      SUBROUTINE OPTPOS
      
      IMPLICIT NONE

      INCLUDE 'cmem.inc'
      INCLUDE 'cout.inc'
      INCLUDE 'cerr.inc'
      
      INCLUDE 'cns.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmg.inc'
      INCLUDE 'cmgpar.inc'
      
      INCLUDE 'cinidat.inc'
      INCLUDE 'ciniopt.inc'
      
      END
            
