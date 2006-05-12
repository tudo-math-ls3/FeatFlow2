***********************************************************************
      SUBROUTINE FPOST(ITYP,IFILEN,ITFILM,UE,MSHOW,DFWX,DFWY,
     *                 ITOPT,COFN)
      
************************************************************************
*   Purpose: - performs the nonsteady postprocess:
*                   - output for film generation
*                   - output for AVS
*                   - output for GMV
*
* In:
*  ITYP   - Type of postprocessing:
*           0=postprocessing of stationary solution
*           1=postprocessing of time-dependent solution
*  IFILEN - Current file number in case of postprocessing of
*           a stationary solution. Is incremented by 1.
*           Not used in case of time-dependent case.
*  ITFILM - Current file number in case of postprocessing of
*           a time-dependent solution. Is incremented by 1.
*           Not used in case of stationary case.
*  UE     - Reference solution for error analysis
*  MSHOW  - Whether to write the results of the analysis to
*           screen/file.
*           0=no output
*           1=write output to file
*           2=write output to file and screen
*  ITOPT  - if <> 0: write current solution into GMV-file with
*                    filename COFN.
*
* Out:
*  DFWX   - Drag-coefficient (if calculated)
*  DFWY   - Lift-coefficient (if calculated)
************************************************************************
      IMPLICIT NONE
      
      INCLUDE 'cmem.inc'
      INCLUDE 'cerr.inc'
      INCLUDE 'cout.inc'
      
      INCLUDE 'cbasicmg.inc'
      INCLUDE 'cmgpar.inc'
      INCLUDE 'cmgfld.inc'
      INCLUDE 'cmg.inc'
      INCLUDE 'cmgtria.inc'
      INCLUDE 'cmgadr.inc'
      
      INCLUDE 'stria.inc'
      INCLUDE 'cmgfbtria.inc'
      
      INCLUDE 'ctria.inc'
      
      INCLUDE 'cns.inc'
      INCLUDE 'cnsparfrac.inc'
      INCLUDE 'cnspts.inc'
      
C *** user COMMON blocks
      INCLUDE 'cinidat.inc'
      INCLUDE 'ciniopt.inc'
      
      INCLUDE 'cfiles.inc'
      
C parameters

      INTEGER ITYP,IFILEN,ITFILM,ITOPT,MSHOW
      CHARACTER*(*) COFN

      DOUBLE PRECISION UE
      EXTERNAL UE

C local variables

      INTEGER NGRAPH,ISETLV,IFMTS,ITWX,KMOV1,KMOV2,KAUXM,LAVSU
      INTEGER KAVS1,KAVS2,IVTA,KPL,LAREA,LAVSP,KAVSP
      INTEGER KISO,KVIND,LAVSI,KAVSI
      INTEGER LAL,LALFA,LMALF,LTMP,I
      
      INTEGER NESTAT(5)

      DOUBLE PRECISION P5, P6, DAW, DFW, DFWX, DFWY, DFWXL, DFWYL
      DOUBLE PRECISION DFRWX, DFRWY
      DOUBLE PRECISION DVOL, DRVOL, DCCDEF(2), DACLNV, DACLNR
      
C=======================================================================
C     Constants
C=======================================================================
      CHARACTER CXX(10)*15,CFILE*15
C
      DATA CXX/'#ns/#DX1       ','#ns/#DX2       ','#ns/#DX3       ',
     *         '#ns/#DX4       ','#ns/#DX5       ','#ns/#DX6       ',
     *         '#ns/#DX7       ','#ns/#DX8       ','#ns/#DX9       ',
     *         '#ns/#DX10      '/

      SAVE CXX

C externals

C Wrapper routine for nonconforming elements. Chooses the correct
C element (E030,E031,EM30,EM31) dependent on the parameter in the
C .DAT file.
      EXTERNAL EXXXNC

      DOUBLE PRECISION FBDVOL
      EXTERNAL FBDVOL
      
C=======================================================================
      SUB='FPOST '

      NGRAPH=MAX(IAVS,IGMV)
      ILEV=NLEV
      ISETLV=2
      CALL  SETLEV (ISETLV)
C
      CALL PTSDAT(TIMENS,DBLE(NY))
C
C=======================================================================
C *** write the solution vector if ABS(ISOL)=1
C=======================================================================
C
      IF ((ITYP.EQ.0).AND.(ABS(ISOL).EQ.1)) THEN
       IF (ISOL.EQ.1) THEN
        IFMTS=0
       ELSE
        IFMTS=1
       ENDIF
       CALL  OF0 (MSOL,CSOL,IFMTS)
       CALL  OWA1 (DWORK(KU1),'DU12P ',NUP,MSOL,IFMTS)
       CLOSE(MSOL)
      ENDIF
C
C=======================================================================
C *** write unformatted time dep. solution vector
C=======================================================================
C
      IF ((ITYP.EQ.1).AND.(INSAV.GT.0)) THEN
      IF (MOD(ITNS,INSAV).EQ.0) THEN
       IFILEN=IFILEN+1
       ITWX=MOD(IFILEN+INSAVN-1,INSAVN)+1
       CALL  OF0 (39,CXX(ITWX),0)
       CALL  OWA1 (DWORK(KU1),'DU12P ',NUP,39,0)
       REWIND(39)
       CLOSE (39)
       IF (IER.NE.0) GOTO 99999
      ENDIF
      ENDIF
C
C=======================================================================
C
      IF ((ITYP.EQ.1).AND.((TIMENS-DTFILO).GE.DTFILM)) THEN
       ITFILM=ITFILM+1
       WRITE(40,*) REAL(TIMENS)
      ENDIF

C=======================================================================

C allocation of temporary vectors for postprocessing

      CALL ZNEW (2*NVT,1,LTMP,'KTMP  ')

C=======================================================================
C *** write velocities 
C=======================================================================

       KMOV1=L(LTMP)
       KMOV2=KMOV1+NVT
       KAUXM=L(KLAUX(NLEV))
C
       CALL  INTUVD(DWORK(KU1),DWORK(KU2),DWORK(KMOV1),DWORK(KMOV2),
     *              DWORK(KAUXM),NVT,NEL,NVBD,KWORK(L(LNPR)),
     *              KWORK(L(LMID)),KWORK(L(LVERT)),KWORK(L(LVBD)),
     *              KWORK(L(KLMBD(NLEV))),DWORK(L(LCORVG)),UE,.FALSE.,
     *              TRIAS(1,ILEV),TIMENS,RE,0,0,0,0)
C     !!!Call incorrect
C
       IF ((ITYP.EQ.1).AND.(IFUSAV.GT.0)
     *                .AND.((TIMENS-DTFILO).GE.DTFILM)) THEN
        CFILE='#film/#DU      '
        IF ((ITFILM.GE.0).AND.(ITFILM.LT.10)) 
     *       WRITE(CFILE(10:10),'(I1.1)') ITFILM
        IF ((ITFILM.GE.10).AND.(ITFILM.LT.100)) 
     *       WRITE(CFILE(10:11),'(I2.2)') ITFILM
        IF ((ITFILM.GE.100).AND.(ITFILM.LT.1000)) 
     *       WRITE(CFILE(10:12),'(I3.3)') ITFILM
        IF ((ITFILM.GE.1000).AND.(ITFILM.LT.10000)) 
     *       WRITE(CFILE(10:13),'(I4.4)') ITFILM
        IF (ITFILM.GE.10000) STOP
C
        CALL OF0 (39,CFILE,0)
        IF (IER.NE.0) GOTO 99998
        CALL  OWA1 (DWORK(KMOV1),'DUF   ',KNVT(IFUSAV),39,0)
        CALL  OWA1 (DWORK(KMOV2),'DUF   ',KNVT(IFUSAV),39,0)
        REWIND(39)
        CLOSE (39)
        IF (IER.NE.0) GOTO 99999
       ENDIF
C
       IF (((ITYP.EQ.0).AND.(IAVS.GT.0)).OR.
     *     ((ITYP.EQ.0).AND.(IGMV.GT.0)).OR.
     *     ((ITYP.EQ.1).AND.(IAVS.GT.0)
     *                 .AND.((TIMENS-DTAVSO).GE.DTAVS)) .OR.
     *     ((ITYP.EQ.1).AND.(IGMV.GT.0)
     *                 .AND.((TIMENS-DTGMVO).GE.DTGMV))) THEN
        CALL ZNEW (2*KNVT(NGRAPH),-2,LAVSU,'AVSU  ')
        IF (IER.NE.0) GOTO 99998
        KAVS1=L(LAVSU)
        KAVS2=KAVS1+KNVT(NGRAPH)
        DO 100 IVTA=1,KNVT(NGRAPH)
        VWORK(KAVS1+IVTA-1)=REAL(DWORK(KMOV1+IVTA-1))
100     VWORK(KAVS2+IVTA-1)=REAL(DWORK(KMOV2+IVTA-1))
       ENDIF
C
      IF (ITYP.EQ.1) THEN
       IF (MSHOW.GE.2) 
     *  WRITE(MTERM,1001) DWORK(KMOV1+KPU(1)-1),DWORK(KMOV2+KPU(1)-1),
     *                    DWORK(KMOV1+KPU(2)-1),DWORK(KMOV2+KPU(2)-1)
       IF (MSHOW.GE.1) 
     *  WRITE(MFILE1,1001) DWORK(KMOV1+KPU(1)-1),DWORK(KMOV2+KPU(1)-1),
     *                    DWORK(KMOV1+KPU(2)-1),DWORK(KMOV2+KPU(2)-1)
       WRITE(41,*) REAL(TIMENS),REAL(DWORK(KMOV1+KPU(1)-1))
       WRITE(42,*) REAL(TIMENS),REAL(DWORK(KMOV2+KPU(1)-1))
       WRITE(43,*) REAL(TIMENS),REAL(DWORK(KMOV1+KPU(2)-1))
       WRITE(44,*) REAL(TIMENS),REAL(DWORK(KMOV2+KPU(2)-1))
      ENDIF
C
C=======================================================================
C *** write PRESSURE 
C=======================================================================

       KPL  =L(LTMP)
       KAUXM=L(LTMP)+NVT
       LAREA=KLAREA(NLEV)
C
       CALL  INTPV (DWORK(KP),DWORK(KPL),DWORK(KAUXM),DWORK(L(LAREA)),
     *              KWORK(L(LVERT)),NVT,NEL)
       IF (IER.NE.0) GOTO 99999
C
       IF ((ITYP.EQ.1).AND.(IFPSAV.GT.0)
     *                .AND.((TIMENS-DTFILO).GE.DTFILM)) THEN
        CFILE='#film/#DP      '
        IF ((ITFILM.GE.0).AND.(ITFILM.LT.10)) 
     *       WRITE(CFILE(10:10),'(I1.1)') ITFILM
        IF ((ITFILM.GE.10).AND.(ITFILM.LT.100)) 
     *       WRITE(CFILE(10:11),'(I2.2)') ITFILM
        IF ((ITFILM.GE.100).AND.(ITFILM.LT.1000)) 
     *       WRITE(CFILE(10:12),'(I3.3)') ITFILM
        IF ((ITFILM.GE.1000).AND.(ITFILM.LT.10000)) 
     *       WRITE(CFILE(10:13),'(I4.4)') ITFILM
        IF (ITFILM.GE.10000) STOP
C
        CALL OF0 (39,CFILE,0)
        IF (IER.NE.0) GOTO 99998
        CALL  OWA1 (DWORK(KPL),'DPL   ',KNVT(IFPSAV),39,0)
        REWIND(39)
        CLOSE (39)
        IF (IER.NE.0) GOTO 99999
       ENDIF
C
       IF (((ITYP.EQ.0).AND.(IAVS.GT.0)).OR.
     *     ((ITYP.EQ.0).AND.(IGMV.GT.0)).OR.
     *     ((ITYP.EQ.1).AND.(IAVS.GT.0)
     *                 .AND.((TIMENS-DTAVSO).GE.DTAVS)) .OR.
     *     ((ITYP.EQ.1).AND.(IGMV.GT.0)
     *                 .AND.((TIMENS-DTGMVO).GE.DTGMV))) THEN
        CALL ZNEW (KNVT(NGRAPH),-2,LAVSP,'AVSP  ')
        IF (IER.NE.0) GOTO 99998
        KAVSP=L(LAVSP)
        DO 200 IVTA=1,KNVT(NGRAPH)
200     VWORK(KAVSP+IVTA-1)=REAL(DWORK(KPL+IVTA-1))
       ENDIF
C
      IF (ITYP.EQ.1) THEN
       IF (MSHOW.GE.2) 
     *  WRITE(MTERM,1002) DWORK(KPL+KPP(1)-1),DWORK(KPL+KPP(2)-1),
     *                    DWORK(KPL+KPP(3)-1),DWORK(KPL+KPP(4)-1)
       IF (MSHOW.GE.1) 
     *  WRITE(MFILE1,1002) DWORK(KPL+KPP(1)-1),DWORK(KPL+KPP(2)-1),
     *                    DWORK(KPL+KPP(3)-1),DWORK(KPL+KPP(4)-1)
       WRITE(45,*) REAL(TIMENS),REAL(DWORK(KPL+KPP(1)-1))
       WRITE(46,*) REAL(TIMENS),REAL(DWORK(KPL+KPP(2)-1))
       WRITE(47,*) REAL(TIMENS),REAL(DWORK(KPL+KPP(3)-1))
       WRITE(48,*) REAL(TIMENS),REAL(DWORK(KPL+KPP(4)-1))
      ENDIF

      CALL  BDPRES(DWORK(KPL),KWORK(L(LVERT)),KWORK(L(LNPR)),
     *             KWORK(L(LVBD)),DWORK(L(LCORVG)),DWORK(L(LVBDP)),
     *             KNVBD(NLEV),P5,P6)

      IF (ITYP.EQ.1) THEN
       IF (MSHOW.GE.2) 
     *  WRITE(MTERM,1003) P5,P6
       IF (MSHOW.GE.1) 
     *  WRITE(MFILE1,1003) P5,P6
       WRITE(49,*) REAL(TIMENS),REAL(P5)
       WRITE(50,*) REAL(TIMENS),REAL(P6)
      ENDIF
C
C=======================================================================
C *** write lift and drag ------<<<<<<<<< featflow original
C=======================================================================
C
C       IF ((IELT.EQ.0).OR.(IELT.EQ.2)) 
C     *  CALL  BDFORC(DWORK(KU1),DWORK(KU2),DWORK(KP),KWORK(L(LVERT)),
C     *               KWORK(L(LMID)),KWORK(L(LVBD)),KWORK(L(LEBD)),
C     *               KWORK(L(LMM)),DWORK(L(LCORVG)),EM31,DFW,DAW)
CC
C       IF ((IELT.EQ.1).OR.(IELT.EQ.3)) 
C     *  CALL  BDFORC(DWORK(KU1),DWORK(KU2),DWORK(KP),KWORK(L(LVERT)),
C     *               KWORK(L(LMID)),KWORK(L(LVBD)),KWORK(L(LEBD)),
C     *               KWORK(L(LMM)),DWORK(L(LCORVG)),EM30,DFW,DAW)
     
      CALL  BDFORC(DWORK(KU1),DWORK(KU2),DWORK(KP),KWORK(L(LVERT)),
     *               KWORK(L(LMID)),KWORK(L(LVBD)),KWORK(L(LEBD)),
     *               KWORK(L(LMM)),DWORK(L(LCORVG)),EXXXNC,DFW,DAW)

      IF (ITYP.EQ.1) THEN
       IF (MSHOW.GE.2) 
     *      WRITE(MTERM,'(A,2(E16.5))') 'I(FORCE) BDR CP             ',
     *                                  DFW,DAW
cc     *  WRITE(MTERM,1004) DFW,DAW
       IF (MSHOW.GE.1) 
     *      WRITE(MFILE1,'(A,2(E16.5))') 'I(FORCE) BDR CP             ',
     *                                  DFW,DAW
cc     *  WRITE(MFILE1,1004) DFW,DAW
       WRITE(51,*) REAL(TIMENS),REAL(DFW),REAL(DAW)
      ENDIF
C
C=======================================================================
C *** write lift and drag by volume integration
C=======================================================================

      CALL ZNEW (NU,1,LAL,'AL   ')

      CALL  BDFVOL(DWORK(KU1),DWORK(KU2),DWORK(KP),
     *    KWORK(L(LVERT)),KWORK(L(LNPR)),
     *    KWORK(L(LMID)),DWORK(L(LCORVG)),EXXXNC,
     *    DWORK(L(lal)),DFWX,DFWY,0)

      CALL ZDISP (0,lal,'al   ')

      IF (ITYP.EQ.1) THEN
       IF (MSHOW.GE.2) 
     *  WRITE(MTERM,'(A,2(E16.5))') 'I(FORCE) VOL CP             ',
     *                              DFWX,DFWY
       IF (MSHOW.GE.1) 
     *  WRITE(MFILE1,'(A,2(E16.5))') 'I(FORCE) VOL CP             ',
     *                              DFWX,DFWY
       WRITE(52,*) REAL(TIMENS),REAL(DFWX),REAL(DFWY)
      ENDIF

C=======================================================================
C *** write lift and drag by volume integration
C More general evaluation method
C=======================================================================

      CALL  BDFVLG(DWORK(KU1),DWORK(KU2),DWORK(KP),DWORK(L(LCORVG)),
     *    KWORK(L(LVERT)),DWORK(L(LCORMG)),KWORK(L(LMID)),
     *    KWORK(L(LNPR)),NEL,NVE,NVT,NBCT,
     *    EXXXNC,(IELT.EQ.2).OR.(IELT.EQ.3),DFWX,DFWY,DPF(1),DPF(2),
     *    NESTAT,0,DVIARX)

      IF (ITYP.EQ.1) THEN
       IF (MSHOW.GE.2) 
     *  WRITE(MTERM,'(A,2(E16.5))') 'I(FORCE) gen. VOL CP        ',
     *                               DFWX,DFWY
       IF (MSHOW.GE.1) 
     *  WRITE(MFILE1,'(A,2(E16.5))') 'I(FORCE) gen. VOL CP        ',
     *                               DFWX,DFWY
       WRITE(52,*) REAL(TIMENS),REAL(DFWX),REAL(DFWY)
      ENDIF
C=======================================================================
C *** write lift and drag by volume integration
C=======================================================================
C
      CALL ZNEW (nu,1,lal,'al   ')

      CALL  BDFVOL(DWORK(KU1),DWORK(KU2),DWORK(KPL),
     *    KWORK(L(LVERT)),KWORK(L(LNPR)),
     *    KWORK(L(LMID)),DWORK(L(LCORVG)),EXXXNC,
     *    DWORK(L(lal)),DFWXL,DFWYL,1)

      CALL ZDISP (0,LAL,'AL   ')

      IF (ITYP.EQ.1) THEN
       IF (MSHOW.GE.2) 
     *  WRITE(MTERM,'(A,2(E16.5))') 'I(FORCE) VOL LP             ',
     *                              DFWXL,DFWYL
       IF (MSHOW.GE.1) 
     *  WRITE(MFILE1,'(A,2(E16.5))') 'I(FORCE) VOL LP             ',
     *                              DFWXL,DFWYL
      ENDIF
      
C=======================================================================
C *** write lift and drag by volume, L2-projected alpha-vector
C     integration
C=======================================================================

       CALL ZNEW(NMT,1,LALFA,'DALFA  ')
       CALL ZNEW(NMT,1,LMALF,'DMALF  ') 

       CALL BDFVL2 (DWORK(KU1),DWORK(KU2),DWORK(KP),NY,
     *                 KWORK(L(LVERT)),KWORK(L(LMID)),KWORK(L(LNPR)),
     *                 DWORK(L(LCORVG)),DWORK(L(LCORMG)),  
     *                 EXXXNC,EXXXNC,IELT,IELT,8,DFWXL,DFWYL,
     *                 DWORK(L(LALFA)),DWORK(L(LMALF)))

       CALL ZDISP(0,LALFA,'DALFA  ')
       CALL ZDISP(0,LMALF,'DMALF  ') 

      IF (ITYP.EQ.1) THEN
       IF (MSHOW.GE.2) 
     *  WRITE(MTERM,'(A,2(E16.5))') 'I(FORCE) VOL LP, L2-ALPHA   ',
     *                              DFWXL,DFWYL
       IF (MSHOW.GE.1) 
     *  WRITE(MFILE1,'(A,2(E16.5))') 'I(FORCE) VOL LP, L2-ALPHA   ',
     *                              DFWXL,DFWYL
      ENDIF

C=======================================================================
C Calculate the vector using the reconstructed line integration method.
C=======================================================================

      CALL  BDFRIG(DWORK(KU1),DWORK(KU2),DWORK(KP),
     *    DWORK(L(LCORVG)),KWORK(L(LVERT)),DWORK(L(LCORMG)),
     *    KWORK(L(LMID)),KWORK(L(LNPR)),NEL,NVE,EXXXNC,
     *    (IELT.EQ.2).OR.(IELT.EQ.3),
     *    DFRWX,DFRWY,DPF(1),DPF(2),0,0,0)
C     Call not correct

      IF (ITYP.EQ.1) THEN
       IF (MSHOW.GE.2) THEN
          WRITE(MTERM,'(A,2(E16.5))') 'I(FORCE) LINE CP, rec. bd.: ',
     *                              DFRWX,DFRWY
       END IF
       IF (MSHOW.GE.1) THEN
         WRITE(MFILE1,'(A,2(E16.5))') 'I(FORCE) LINE CP, rec. bd.: ',
     *                               DFRWX,DFRWY
       END IF
      ENDIF

C=======================================================================
C calculation of the volume that is covered by fictitious boundary
C components
C=======================================================================
      
      CALL AFBVOL (KWORK(L(LVERT)),KWORK(L(LMID)),
     *               DWORK(L(LCORVG)),EXXXNC,8,0,DVOL)

C=======================================================================
C Calculation of the Closed-Curve-Defect
C
C This performs an integration using the prescribed vector
C (U1,U2,P)=(0,0,1), so the surface integral about the circle
C (to be more exact: about a closed curve of the boundary of an
C object) 
C        int_S ( [grad(U) + p*I] * normal_vec ) ds
C which is used for the calculation of drag/lift on the surface of
C the obstacle in the channel becomes the zero vector - theoretically!
C The difference-vector to zero is called "Closed-Curve-Defect" here
C and indicates how accurate the boundary is approximated by the
C grid. This vector should tend to zero with further refinements,
C expecially if adaptive grid refinement is used!
C=======================================================================

C allocate auxiliary "solution" vector
      
      CALL ZNEW(2*NMT+NEL,1,LAL,'DAUSOL')
      
C fill it with data

      DO I=2*NMT,2*NMT+NEL-1
        DWORK(L(LAL)+I)=1
      END DO

C Calculate the vector using the general volume integration method.
C There is no drag/lift correction involved here, so DPF1=DPF2=1.

      CALL  BDFVLG(DWORK(L(LAL)),DWORK(L(LAL)+NMT),DWORK(L(LAL)+2*NMT),
     *    DWORK(L(LCORVG)),KWORK(L(LVERT)),DWORK(L(LCORMG)),
     *    KWORK(L(LMID)),KWORK(L(LNPR)),NEL,NVE,NVT,NBCT,
     *    EXXXNC,(IELT.EQ.2).OR.(IELT.EQ.3),
     *    DCCDEF(1),DCCDEF(2),1D0,1D0,NESTAT,0,DVIARX)

C Clean up memory

      CALL ZDISP (0,LAL,'DAUSOL')

C=======================================================================
C Calculation of approximate arc length of fictitious boundary
C
C This uses a special trick to calculate the arc length of all
C fictitious boundary components with a volume integration approach,
C see ARCLENGTH.F.
C=======================================================================

C Calculate the vector using the general volume integration method.
C There is no drag/lift correction involved here, so DPF1=DPF2=1.

      DACLNV=0D0
      DACLNR=0D0

      CALL  BDARCL(DWORK(L(LCORVG)),KWORK(L(LVERT)),
     *          DWORK(L(LCORMG)),KWORK(L(LMID)),KWORK(L(LNPR)),
     *          NEL,NVE,NVT,NBCT,EXXXNC,0,DACLNV, DACLNR,DVIARX,0,0)
      !!! call not correct

C=======================================================================
C print the results
C=======================================================================

C reference volume

      DRVOL = FBDVOL (0)

      IF (ITYP.EQ.1) THEN
        IF (MSHOW.GE.2) THEN
          WRITE(MTERM,'(A,(D24.14))') 'approx. VOL(FBDRY) : ',DVOL
          WRITE(MTERM,'(A,(D20.10))') 'ref. VOL(FBDRY)    : ',DRVOL
          WRITE(MTERM,'(A,(D20.10))') 'VOL-eff.(FBDRY)    : ',DVOL/DRVOL
          WRITE(MTERM,'(A,(2D20.10))') 
     *         'Closed-Curve-Defect:     ',DCCDEF(1),DCCDEF(2)
          WRITE(MTERM,'(A,(2D20.10))') 
     *         'ARC-Length of fict. bd.: ',DACLNV, DACLNR
          WRITE(MTERM,'(A,I7,A,D18.8)') 
     *      '#Elements with 0 edge-midpoint in fict. bd.: ',NESTAT(1),
     *      ' /NEL = ',DBLE(NESTAT(1))/DBLE(NEL)
          WRITE(MTERM,'(A,I7,A,D18.8)') 
     *      '#Elements with 1 edge-midpoint in fict. bd.: ',NESTAT(2),
     *      ' /NEL = ',DBLE(NESTAT(2))/DBLE(NEL)
          WRITE(MTERM,'(A,I7,A,D18.8)') 
     *      '#Elements with 2 edge-midpoint in fict. bd.: ',NESTAT(3),
     *      ' /NEL = ',DBLE(NESTAT(3))/DBLE(NEL)
          WRITE(MTERM,'(A,I7,A,D18.8)') 
     *      '#Elements with 3 edge-midpoint in fict. bd.: ',NESTAT(4),
     *      ' /NEL = ',DBLE(NESTAT(4))/DBLE(NEL)
          WRITE(MTERM,'(A,I7,A,D18.8)') 
     *      '#Elements with 4 edge-midpoint in fict. bd.: ',NESTAT(5),
     *      ' /NEL = ',DBLE(NESTAT(5))/DBLE(NEL)
        END IF
        IF (MSHOW.GE.1) THEN
          WRITE(MFILE1,'(A,(D20.10))') 'approx. VOL(FBDRY) : ',DVOL
          WRITE(MFILE1,'(A,(D20.10))') 'ref. VOL(FBDRY)    : ',DRVOL
          WRITE(MFILE1,'(A,(D20.10))') 'VOL-eff.(FBDRY)    : ',
     *                                 DVOL/DRVOL
          WRITE(MFILE1,'(A,(2D20.10))') 
     *         'Closed-Curve-Defect: ',DCCDEF(1),DCCDEF(2)
          WRITE(MFILE1,'(A,(2D20.10))') 
     *         'ARC-Length of fict. bd.: ',DACLNV, DACLNR
          WRITE(MFILE1,'(A,I7,A,D18.8)') 
     *      '#Elements with 0 edge-midpoint in fict. bd.: ',NESTAT(1),
     *      ' /NEL = ',DBLE(NESTAT(1))/DBLE(NEL)
          WRITE(MFILE1,'(A,I7,A,D18.8)') 
     *      '#Elements with 1 edge-midpoint in fict. bd.: ',NESTAT(2),
     *      ' /NEL = ',DBLE(NESTAT(2))/DBLE(NEL)
          WRITE(MFILE1,'(A,I7,A,D18.8)') 
     *      '#Elements with 2 edge-midpoint in fict. bd.: ',NESTAT(3),
     *      ' /NEL = ',DBLE(NESTAT(3))/DBLE(NEL)
          WRITE(MFILE1,'(A,I7,A,D18.8)') 
     *      '#Elements with 3 edge-midpoint in fict. bd.: ',NESTAT(4),
     *      ' /NEL = ',DBLE(NESTAT(4))/DBLE(NEL)
          WRITE(MFILE1,'(A,I7,A,D18.8)') 
     *      '#Elements with 4 edge-midpoint in fict. bd.: ',NESTAT(5),
     *      ' /NEL = ',DBLE(NESTAT(5))/DBLE(NEL)
        END IF
      ENDIF
      
C=======================================================================
C *** write streamfunction 
C=======================================================================

       KISO =L(LTMP)
       KVIND=L(LTMP)+NVT
     
       CALL  U2ISO (DWORK(L(LCORVG)),KWORK(L(LVERT)),KWORK(L(LMID)),
     *              KWORK(L(LADJ)),NVT,NEL,NVE,
     *              DWORK(KVIND),DWORK(KISO),
     *              DWORK(KU1),DWORK(KU2))

       IF (IER.NE.0) GOTO 99999
C
       IF ((ITYP.EQ.1).AND.(IFXSAV.GT.0)
     *                .AND.((TIMENS-DTFILO).GE.DTFILM)) THEN
        CFILE='#film/#DX      '
        IF ((ITFILM.GE.0).AND.(ITFILM.LT.10)) 
     *       WRITE(CFILE(10:10),'(I1.1)') ITFILM
        IF ((ITFILM.GE.10).AND.(ITFILM.LT.100)) 
     *       WRITE(CFILE(10:11),'(I2.2)') ITFILM
        IF ((ITFILM.GE.100).AND.(ITFILM.LT.1000)) 
     *       WRITE(CFILE(10:12),'(I3.3)') ITFILM
        IF ((ITFILM.GE.1000).AND.(ITFILM.LT.10000)) 
     *       WRITE(CFILE(10:13),'(I4.4)') ITFILM
        IF (ITFILM.GE.10000) STOP
C
        CALL OF0 (39,CFILE,0)
        IF (IER.NE.0) GOTO 99998
        CALL  OWA1 (DWORK(KISO),'DISO  ',KNVT(IFXSAV),39,0)
        REWIND(39)
        CLOSE (39)
        IF (IER.NE.0) GOTO 99999
       ENDIF
C
       IF (((ITYP.EQ.0).AND.(IAVS.GT.0)).OR.
     *     ((ITYP.EQ.0).AND.(IGMV.GT.0)).OR.
     *     ((ITYP.EQ.1).AND.(IAVS.GT.0)
     *                 .AND.((TIMENS-DTAVSO).GE.DTAVS)) .OR.
     *     ((ITYP.EQ.1).AND.(IGMV.GT.0)
     *                 .AND.((TIMENS-DTGMVO).GE.DTGMV))) THEN
        CALL ZNEW (KNVT(NGRAPH),-2,LAVSI,'AVSI  ')
        IF (IER.NE.0) GOTO 99998
        KAVSI=L(LAVSI)
        DO 300 IVTA=1,KNVT(NGRAPH)
300     VWORK(KAVSI+IVTA-1)=REAL(DWORK(KISO+IVTA-1))
       ENDIF
C
      IF (ITYP.EQ.1) THEN
       IF (MSHOW.GE.2) 
     *  WRITE(MTERM,1005) DWORK(KISO+KPX(1)-1)-DWORK(KISO+KPX(2)-1)
       IF (MSHOW.GE.1) 
     *  WRITE(MFILE1,1005) DWORK(KISO+KPX(1)-1)-DWORK(KISO+KPX(2)-1)
       WRITE(53,*) REAL(TIMENS),REAL( DWORK(KISO+KPX(1)-1)
     *                               -DWORK(KISO+KPX(2)-1))
      ENDIF
C
C=======================================================================
C
      IF ((ITYP.EQ.1).AND.((TIMENS-DTFILO).GE.DTFILM)) THEN
       DTFILO=TIMENS
      ENDIF
C
C=======================================================================
C *** write AVS 
C=======================================================================
C
      IF (((ITYP.EQ.0).AND.(IAVS.GT.0)).OR.
     *    ((ITYP.EQ.0).AND.(IGMV.GT.0)).OR.
     *    ((ITYP.EQ.1).AND.(IAVS.GT.0)
     *                .AND.((TIMENS-DTAVSO).GE.DTAVS)) .OR.
     *    ((ITYP.EQ.1).AND.(IGMV.GT.0)
     *                .AND.((TIMENS-DTGMVO).GE.DTGMV))) THEN
C
       IF (((ITYP.EQ.0).AND.(IAVS.GT.0)).OR.
     *     ((ITYP.EQ.1).AND.(IAVS.GT.0)
     *                 .AND.((TIMENS-DTAVSO).GE.DTAVS))) THEN
       CFILE='#avs/u.        '
       IF ((ITNS+IFINIT.GE.0).AND.(ITNS+IFINIT.LT.10)) 
     *      WRITE(CFILE(8:8),'(I1.1)') ITNS+IFINIT
       IF ((ITNS+IFINIT.GE.10).AND.(ITNS+IFINIT.LT.100)) 
     *      WRITE(CFILE(8:9),'(I2.2)') ITNS+IFINIT
       IF ((ITNS+IFINIT.GE.100).AND.(ITNS+IFINIT.LT.1000)) 
     *      WRITE(CFILE(8:10),'(I3.3)') ITNS+IFINIT
       IF ((ITNS+IFINIT.GE.1000).AND.(ITNS+IFINIT.LT.10000)) 
     *      WRITE(CFILE(8:11),'(I4.4)') ITNS+IFINIT
       IF (ITNS+IFINIT.GE.10000) STOP

       CALL XAVS2D(39,CFILE,KNEL(IAVS),KNVT(IAVS),
     *             KWORK(L(KLVERT(IAVS))),DWORK(L(LCORVG)),
     *             VWORK(KAVS1),VWORK(KAVS2),VWORK(KAVSP),
     *             VWORK(KAVSI))
       DTAVSO=TIMENS
       ENDIF
C
C=======================================================================
C *** write GMV 
C=======================================================================

       IF (((ITYP.EQ.0).AND.(IGMV.GT.0)).OR.
     *     ((ITYP.EQ.1).AND.(IGMV.GT.0)
     *                 .AND.((TIMENS-DTGMVO).GE.DTGMV))) THEN
         CFILE='#gmv/u.        '       
         IF ((ITNS+IFINIT.GE.0).AND.(ITNS+IFINIT.LT.10)) 
     *      WRITE(CFILE(8:8),'(I1.1)') ITNS+IFINIT
         IF ((ITNS+IFINIT.GE.10).AND.(ITNS+IFINIT.LT.100)) 
     *      WRITE(CFILE(8:9),'(I2.2)') ITNS+IFINIT
         IF ((ITNS+IFINIT.GE.100).AND.(ITNS+IFINIT.LT.1000)) 
     *      WRITE(CFILE(8:10),'(I3.3)') ITNS+IFINIT
         IF ((ITNS+IFINIT.GE.1000).AND.(ITNS+IFINIT.LT.10000)) 
     *      WRITE(CFILE(8:11),'(I4.4)') ITNS+IFINIT
         IF (ITNS+IFINIT.GE.10000) STOP
 
         CALL XGMV2D(39,CFILE,KNEL(IGMV),KNVT(IGMV),
     *             KWORK(L(KLVERT(IGMV))),DWORK(L(LCORVG)),
     *             KWORK(L(KLNPR(IGMV))),
     *             VWORK(KAVS1),VWORK(KAVS2),VWORK(KAVSP),VWORK(KAVSI),
     *             DWORK(L(KLAREA(IGMV))),TIMENS)


         IF (ITOPT.NE.0) THEN

           CALL XGMV2D(39,COFN,KNEL(IGMV),KNVT(IGMV),
     *             KWORK(L(KLVERT(IGMV))),DWORK(L(LCORVG)),
     *             KWORK(L(KLNPR(IGMV))),
     *             VWORK(KAVS1),VWORK(KAVS2),VWORK(KAVSP),VWORK(KAVSI),
     *             DWORK(L(KLAREA(IGMV))),TIMENS)
         END IF

         DTGMVO=TIMENS
       END IF

       CALL ZDISP(0,LAVSI,'AVSI  ')
       CALL ZDISP(0,LAVSP,'AVSP  ')
       CALL ZDISP(0,LAVSU,'AVSU  ')
       IF (IER.NE.0) GOTO 99998
      END IF
C
C=======================================================================
C Release temporary vector

      CALL ZDISP (0,LTMP,'KTMP  ')

      GOTO 99999
C
C=======================================================================
C     Error case
C=======================================================================
99998 WRITE(MTERM,*) 'IER', IER
      WRITE(MTERM,*) 'IN SUBROUTINE ',SUB
C
1000  FORMAT (6E12.5)
1001  FORMAT ('P(VELO) ',4(D12.5))
1002  FORMAT ('P(PRES) ',4(D12.5))
1102  FORMAT ('P(DIFF) ',1(D12.5))
1003  FORMAT ('I(PRES) ',2(D12.5))
1004  FORMAT ('I(FORCE)',2(D12.5))
1005  FORMAT ('P(FLUX) ',1(D12.5))
c
c
c
99999 END
