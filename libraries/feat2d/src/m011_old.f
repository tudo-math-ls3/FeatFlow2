************************************************************************
*    M011   (edited from M010)                                         *
************************************************************************
* FINITE ELEMENT ANALYSIS TOOLBOX  FEAT  (Release 1.3)                 *
*                                                                      *
* Authors: H. Blum, J. Harig, S. Mueller, S. Turek                     *
*          Institute of Applied Mathematics                            *
*          University of Heidelberg                                    *
*          D-6900 HEIDELBERG                                           *
*                                                                      *
************************************************************************
*                                                                      *
* M011                                                                 *
*                                                                      *
* Purpose  Solution of a linear system  A*X = B  using                 *
*          multigrid iteration                                         *
*          Double precision version                                    *
*                                                                      *
* Subroutines/functions called  LL21 , LLC1                            *
*                                                                      *
* Version from  08/25/90                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* DX       R*8    Starting address of vectors containing the           *
* DB       R*8    solution and the right hand side, DD is used as      *
* DD       R*8    auxiliary vector only                                *
* KOFFX           The actual starting address of DX on level ILEV      *
* KOFFB           is DX(1+KOFFX(ILEV)) (analogously for DB and DD)     *
* KOFFD           Total space required for all vectors is              *
*                 KNEQ(NLMIN)+...+KNEQ(NLMAX)                          *
*                 DX(1+KOFFX(NLMAX)) contains initial solution         *
*                 DB(1+KOFFB(NLMAX)) contains right hand side          *
* KNEQ     I*4    Number of equations for all levels                   *
* NLMAX    I*4    Iteration uses levels NLMIN to NLMAX,                *
* NLMIN    I*4    NLMAX  is the finest level                           *
* NIT      I*4    Maximum number of iterations                         *
*                 Iteration completed after reaching the finest level  *
* EPS1      R*8    Desired precision                                   *
*                 Stop if !!DEFN!! < EPS1 !!DEF0!!                    *
* EPS2      R*8    Desired precision                                   *
*                 Stop if !!DEF!! < EPS2                               *
* KPRSM    I*4    Number of pre -smoothing steps for all levels        *
* KPOSM    I*4    Number of post-smoothing steps for all levels        *
* ICYCLE   I*4    <0: special cycle types (not yet implemented)        *
*                 =0: F-Cycle                                          *
*                 =1: V-Cycle                                          *
*                 =2: W-Cycle                                          *
*                 >2: Cycle of higher order                            *
* DAX      SUBR   CALL DAX(DX,DAX,NEQ,A1,A2)                           *
*                 Returns  DAX := A1*A*DX+A2*DAX                       *
* DPROL    SUBR   DPROL(DX1,DX2)                                       *
*                 Returns  DX2 := Prolongation(DX1)  to higher level   *
* DSTEP    SUBR   DPROL(DX,DD,DB,DSTEPP)                               *
*                 Returns  DSTEPP := optimal step size for correction  *
* DREST    SUBR   DREST(DD1,DD2)                                       *
*                 Returns  DD2 := Restriction(DD1)  to lower level     *
* DPRSM    SUBR   DPRSM(DX,DB,DD,NEQ,NPRSM)                            *
*                 Returns  DX  after  NPRSM:=KPRSM(ILEV)               *
*                 pre-smoothing steps, DD  is used as auxiliary vector *
* DPOSM    SUBR   Same as above, used for post-smoothing               *
* DEX      SUBR   DEX(DX,DB,DD,NEQ)                                    *
*                 Returns exact solution                               *
* DBC      SUBR   DBC(DX,NEQ)                                          *
*                 Copies boundary data onto components of  DX          *
* KIT0     I*4    auxiliary vectors of length  NLMAX                   *
* KIT      I*4                                                         *
* BNRM     BOL    Check residuum in euclidian norm. If set to .FALSE., *
*                 the norm is measured in the (weighted) l2-norm.      *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* DX       R*8    Solution vector on DX(KOFFX(NLMAX))                  *
* ITE      I*4    Number of iterations                                 *
* IER      I*4    Error indicator                                      *
* DEF      R*8    final defect                                         *
* RHOLMG   R*8    Multigrid convergence rate                           *
* RHOASM   R*8    MULTIGRID ASYMPTOTIC CONVERGENCE RATE                *
*                 calculated from the last three residuals             *
* BMGEND   BOL    is set to .TRUE. if MG did not converge / error      *
*                                                                      *
************************************************************************
C
      SUBROUTINE  M011 (DX,DB,DD,KOFFX,KOFFB,KOFFD,KNEQ,NIT,ITE,
     *                  EPS1,EPS2,DEF,DAX,DPROL,DREST,DPRSM,DPOSM,
     *                  DEX,DEXA,DBC,DSTEP,KIT0,KIT,IREL,IDEFMG,RHOLMG,
     *                  RHOASM,BMGEND,BNRM)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
C
      PARAMETER (NNARR=299,NNLEV=9)
      
C A possibility to skip a level in multigrid:
C Enter here a level that should be skipped in prolongation/restriction.
C If set to -1, no level is skipped.
      PARAMETER (ISKPLV=-1)
      
      DIMENSION DX(*),DB(*),DD(*),KOFFX(*),KOFFB(*),KOFFD(*)
      DIMENSION KNEQ(*),KIT0(*),KIT(*)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /OUTPUT/ M,MT,MKEYB,MTERM,MERR,MPROT,MSYS,MTRC,IRECL8
      COMMON /MGTRD/  KNEL(NNLEV),KNVT(NNLEV),KNMT(NNLEV),
     *                KNVEL(NNLEV),KNVBD(NNLEV)
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGTIME/ TTMG,TTS,TTE,TTD,TTP,TTR,IMTIME
      EXTERNAL DAX,DPROL,DREST,DPRSM,DPOSM,DEX,DEXA,DBC,DSTEP

C Debug:
C      PARAMETER (NNARR=299)
C      COMMON          NWORK,IWORK,IWMAX,L(299),DWORK(1)
C      COMMON /MGFLD/  KLA(NNLEV),KLCOLA(NNLEV),
C     *                KLLDA(NNLEV),KLU(NNLEV),
C     *                KLB(NNLEV),KLD(NNLEV),KLAUX(NNLEV)
C      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C Debug Ende

      SAVE /ERRCTL/,/CHAR/,/OUTPUT/,/MGTRD/,/MGPAR/,/MGTIME/
      
C Length of the queue of last residuals for the computation of
C the asymptotic convergence rate
      PARAMETER (IASRLN=3)
C The queue saves the current residual and the two previous
C and the two previous residuals
      DIMENSION RESQUE(IASRLN)
      DOUBLE PRECISION RESQUE

C Debug: Ausgabe der Anfangsdaten auf feinstem Level
C      print *,'Startvektor:'
C      write (*,'(80D26.18)') (dx(koffx(NLMAX)+i),i=1,80)
C      print *,'Rechte Seite:'
C      write (*,'(80D26.18)') (db(koffb(NLMAX)+i),i=1,80)
C      print *,'Matrix:'
C      write (*,'(80D26.18)') (dwork(L(KLA(NLMAX))+i-1),i=1,80)
C Debug Ende

      SUB='M011  '
      IF (ICHECK.GE.997) CALL OTRC('M011  ','05/13/91')
      IER=0
C
      BMGEND=.FALSE.
      BREL=IREL.EQ.1
      BMSG2=M.GE.2.OR.MT.GE.2
      MT0=MT
      MT=0
C
      BTIME=IMTIME.GT.0
      IF (BTIME) THEN
       IF (IMTIME.EQ.1) THEN
        TTMG=0D0
        TTS=0D0
        TTE=0D0
        TTD=0D0
        TTP=0D0
        TTR=0D0
       ENDIF
       CALL ZTIME(TTMG0)
      ENDIF
C
C
      NIT0=MAX(ITE,0)
      IF (ICYCLE.LT.0) THEN
       CALL WERR(-181,'M011  ')
       GOTO 99999
      ENDIF
      IF (NLMAX.LT.NLMIN.OR.NLMIN.LE.0.OR.NLMAX.GT.NNLEV) THEN
       CALL WERR(-182,'M011  ')
       GOTO 99999
      ENDIF
      ILEV=NLMAX
C
C     WIR FANGEN MIT NULL AN!
C
CCCCC      CALL LCL1(DX(1+KOFFX(NLMAX)),KNEQ(NLMAX))
C
C
C
C *** special case - zero solution
      IF (BTIME) CALL ZTIME(TTD0)
      CALL LLI1(DB(1+KOFFB(NLMAX)),KNEQ(NLMAX),DMAX,INDMAX)
      IF (DMAX.LE.1D-12) THEN
       CALL LCP1(DB(1+KOFFB(NLMAX)),DX(1+KOFFX(NLMAX)),KNEQ(NLMAX))
       RHOLMG=0D0
       GOTO 99999
      ENDIF
C
      IF (BTIME) THEN
       CALL ZTIME(TTD1)
       TTD=TTD+TTD1-TTD0
      ENDIF
C
C *** special case - only one level
      IF (NLMIN.EQ.NLMAX) THEN
       IF (MT0.GE.2) THEN
         WRITE (MTERM,'(A)') ' *   M011     *   Only one level;'//
     *                       ' switching back to standard solver.'
       END IF
       IF (BTIME) CALL ZTIME(TTE0)
       CALL DEXA(DX(1+KOFFX(NLMAX)),DB(1+KOFFB(NLMAX)),
     *           DD(1+KOFFD(NLMAX)),KNEQ(NLMAX),RHOLMG,EPS2,ITECG)
       CALL LCP1(DB(1+KOFFB(NLMAX)),DD(1+KOFFD(NLMAX)),KNEQ(NLMAX))
       CALL DAX (DX(1+KOFFX(NLMAX)),DD(1+KOFFD(NLMAX)),KNEQ(NLMAX),
     *           -1D0,1D0)
       CALL LL21(DD(1+KOFFD(NLMAX)),KNEQ(NLMAX),DEF)
C Scaling for the vector (1111...) to have norm 1 (weighted l2-norm)
       IF (.NOT.BNRM) DEF = DEF / SQRT (DBLE(KNEQ(NLMAX)))
       IF (BTIME) THEN
        CALL ZTIME(TTE1)
        TTE=TTE+TTE1-TTE0
       ENDIF

C In this case we can't get the asymptotic residual, otherwise we would
C have to change the code of DEXA.
       RHOASM = RHOLMG
C initialize the queue of the last residuals with DEF
       DO I=1,IASRLN
         RESQUE(I)=DEF
       END DO  

       ITE=ITECG
       GOTO 99999
      ENDIF
C
C *** level counts
      KIT0(NLMAX)=1
      DO 2 ILEV=NLMIN+1,NLMAX-1
      IF (ICYCLE.EQ.0) THEN
       KIT0(ILEV)=2
      ELSE
       KIT0(ILEV)=ICYCLE
      ENDIF
2     CONTINUE
C
CCC      IF (BTIME) CALL ZTIME(TTS0)
CCC      IF (KPRSM(NLMAX).GT.0) 
CCC     * CALL DPRSM(DX(1+KOFFX(NLMAX)),DB(1+KOFFB(NLMAX)),
CCC     *            DD(1+KOFFD(NLMAX)),KNEQ(NLMAX),KPRSM(NLMAX))
CCC
CCC      IF (BTIME) THEN
CCC       CALL ZTIME(TTS1)
CCC       TTS=TTS+TTS1-TTS0
CCC      ENDIF
C
      IF (IDEFMG.EQ.1) THEN
       IF (BTIME) CALL ZTIME(TTD0)
       CALL LCP1(DB(1+KOFFB(NLMAX)),DD(1+KOFFD(NLMAX)),KNEQ(NLMAX))
       CALL DAX(DX(1+KOFFX(NLMAX)),DD(1+KOFFD(NLMAX)),KNEQ(NLMAX),
     *          -1D0,1D0)
C
       CALL LL21(DD(1+KOFFD(NLMAX)),KNEQ(NLMAX),DEF)
C Scaling for the vector (1111...) to have norm 1 (weighted l2-norm)
       IF (.NOT.BNRM) DEF = DEF / SQRT (DBLE(KNEQ(NLMAX)))
       DEFOLD=DEF
C *** FD  is considered as initial defect
       FD=DEF

C Initialize the queue of the last residuals with DEF
       DO I=1,IASRLN
         RESQUE(I)=DEF
       END DO  
C
       IF (BTIME) THEN
        CALL ZTIME(TTD1)
        TTD=TTD+TTD1-TTD0
       ENDIF
C
       IF (DEF.LE.EPS2.AND..NOT.BREL) THEN
        ITE=0
        GOTO 1000
       ENDIF
      ENDIF
C
C
C *** Start multigrid iteration
      DO 100 ITE=1,NIT
C
C *** initialize level counts for all levels
      DO 101 ILEV=NLMIN,NLMAX
101   KIT(ILEV)=KIT0(ILEV)
C
      ILEV=NLMAX
110   IF (ILEV.NE.NLMIN) THEN
C
CCC      IF (ILEV.NE.NLMAX) THEN
CCC ***  defect on finest level already available
C
        IF (BTIME) CALL ZTIME(TTS0)
        IF (KPRSM(ILEV).GT.0) 
     *   CALL DPRSM(DX(1+KOFFX(ILEV)),DB(1+KOFFB(ILEV)),
     *              DD(1+KOFFD(ILEV)),KNEQ(ILEV),KPRSM(ILEV))
C
        IF (BTIME) THEN
         CALL ZTIME(TTS1)
         TTS=TTS+TTS1-TTS0
        ENDIF
C
        IF (BTIME) CALL ZTIME(TTD0)
        CALL LCP1(DB(1+KOFFB(ILEV)),DD(1+KOFFD(ILEV)),KNEQ(ILEV))
        CALL DAX(DX(1+KOFFX(ILEV)),DD(1+KOFFD(ILEV)),KNEQ(ILEV),
     *           -1D0,1D0)
C
        IF (BTIME) THEN
         CALL ZTIME(TTD1)
         TTD=TTD+TTD1-TTD0
        ENDIF
C
CCC       ENDIF
C
       ILEV=ILEV-1
C
C *** restriction of defect
       IF (BTIME) CALL ZTIME(TTR0)
       CALL DREST(DD(1+KOFFD(ILEV+1)),DB(1+KOFFB(ILEV)))
       
C Skip Level ISKPLV
        IF (ILEV.EQ.ISKPLV) THEN
          CALL LCP1(DB(1+KOFFB(ILEV)),DD(1+KOFFD(ILEV)),KNEQ(ILEV))
          ILEV = ILEV -1
          CALL DREST(DD(1+KOFFD(ILEV+1)),DB(1+KOFFB(ILEV)))
        END IF
C
       IF (BTIME) THEN
        CALL ZTIME(TTR1)
        TTR=TTR+TTR1-TTR0
       ENDIF
C
C ***  choose zero as initial vector on lower level
       IF (BTIME) CALL ZTIME(TTD0)
       CALL LCL1(DX(1+KOFFX(ILEV)),KNEQ(ILEV))
       CALL DBC(DB(1+KOFFB(ILEV)),KNEQ(ILEV))
       IF (BTIME) THEN
        CALL ZTIME(TTD1)
        TTD=TTD+TTD1-TTD0
       ENDIF
C
       GOTO 110
C
      ENDIF
C
C *** exact solution on lowest level
      IF (BTIME) CALL ZTIME(TTE0)
      CALL DEX(DX(1+KOFFX(NLMIN)),DB(1+KOFFB(NLMIN)),
     *        DD(1+KOFFD(NLMIN)),KNEQ(NLMIN),RHOLMG)
C
      IF (BTIME) THEN
       CALL ZTIME(TTE1)
       TTE=TTE+TTE1-TTE0
      ENDIF
C
C
130   IF (ILEV.NE.NLMAX) THEN
       ILEV=ILEV+1
C
C *** DPROL  returns  DD:=PROL(DX)
      IF (BTIME) CALL ZTIME(TTP0)
      CALL DPROL(DX(1+KOFFX(ILEV-1)),DD(1+KOFFD(ILEV)))

C Skip Level ISKPLV
      IF (ILEV.EQ.ISKPLV) THEN
        CALL LCP1 (DD(1+KOFFD(ILEV)),DX(1+KOFFX(ILEV)),KNEQ(ILEV))
        ILEV = ILEV+1
        CALL DPROL(DX(1+KOFFX(ILEV-1)),DD(1+KOFFD(ILEV)))
      END IF
      
      
      CALL DBC(DD(1+KOFFD(ILEV)),KNEQ(ILEV))
C
      CALL DSTEP(DX(1+KOFFX(ILEV)),DD(1+KOFFD(ILEV)),
     *           DB(1+KOFFB(ILEV)),KNEQ(ILEV),DSTEPP)
      CALL LLC1(DD(1+KOFFD(ILEV)),DX(1+KOFFX(ILEV)),KNEQ(ILEV),
     *          DSTEPP,1D0)
C
      IF (BTIME) THEN
       CALL ZTIME(TTP1)
       TTP=TTP+TTP1-TTP0
      ENDIF
C
C
C *** Post-smoothing
      IF (BTIME) CALL ZTIME(TTS0)
      IF (KPOSM(ILEV).GT.0) 
     * CALL DPOSM(DX(1+KOFFX(ILEV)),DB(1+KOFFB(ILEV)),
     *            DD(1+KOFFD(ILEV)),KNEQ(ILEV),KPOSM(ILEV))
C
      IF (BTIME) THEN
       CALL ZTIME(TTS1)
       TTS=TTS+TTS1-TTS0
      ENDIF
C
       KIT(ILEV)=KIT(ILEV)-1
       IF (KIT(ILEV).EQ.0) THEN
        IF (ICYCLE.EQ.0) THEN
         KIT(ILEV)=1
        ELSE
         KIT(ILEV)=KIT0(ILEV)
        ENDIF
        GOTO 130
       ELSE
        GOTO 110
       ENDIF
C
      ENDIF
C
C
CCC      IF (BTIME) CALL ZTIME(TTS0)
CCC      IF (KPRSM(NLMAX).GT.0)  
CCC     * CALL DPRSM(DX(1+KOFFX(NLMAX)),DB(1+KOFFB(NLMAX)),
CCC     *           DD(1+KOFFD(NLMAX)),KNEQ(NLMAX),KPRSM(NLMAX))
CCC
CCC      IF (BTIME) THEN
CCC       CALL ZTIME(TTS1)
CCC       TTS=TTS+TTS1-TTS0
CCC      ENDIF
C
C
      IF (IDEFMG.EQ.1) THEN
       IF (BTIME) CALL ZTIME(TTD0)
       CALL LCP1(DB(1+KOFFB(NLMAX)),DD(1+KOFFD(NLMAX)),KNEQ(NLMAX))
       CALL DAX(DX(1+KOFFX(NLMAX)),DD(1+KOFFD(NLMAX)),KNEQ(NLMAX),
     *          -1D0,1D0)
C
       CALL LL21(DD(1+KOFFD(NLMAX)),KNEQ(NLMAX),DEF)
C Scaling for the vector (1111...) to have norm 1 (weighted l2-norm)
       IF (.NOT.BNRM) DEF = DEF / SQRT (DBLE(KNEQ(NLMAX)))
       
C shift the queue with the last residuals and add the new
C residual to it
       DO I=1,IASRLN-1
         RESQUE(I)=RESQUE(I+1)
       END DO  
       RESQUE(IASRLN)=DEF     
       
       IF (MT.GE.9) WRITE (6,10001) ITE,DEF
C       IF (MT.GT.0) WRITE (6,10001) ITE,DEF
       IF (BMSG2) THEN
        MT=MT0
        WRITE (CPARAM,'(I15,D25.16)') ITE,DEF
C Michael Koester: If the number is near infinity, don't print it;
C will result in an error
        IF (DEF.LT.1D99) CALL OMSG(73,'M011  ')
C        WRITE (MPROT,10001) ITE,DEF
C        WRITE (*,10001) ITE,DEF
        MT=0
       ENDIF
C
       IF (BTIME) THEN
        CALL ZTIME(TTD1)
        TTD=TTD+TTD1-TTD0
       ENDIF
C
C=======================================================================
C ***  Unexpected STOP !!!
C=======================================================================
       IF ((DEF.GE.1D99).OR.(DEF/DEFOLD.GT.1D3).OR.
     *    (DEF/FD.GT.5D3)) THEN
        BMGEND=.TRUE.
        GOTO 1000
       ENDIF
C
      DEFOLD=DEF
       IF (BREL) THEN
        IF (DEF.LE.FD*EPS1.AND.DEF.LE.EPS2.AND.ITE.GE.NIT0) GOTO 1000
       ELSE
        IF (DEF.LE.EPS2) GOTO 1000
       ENDIF
      ENDIF
C
100   CONTINUE
C
      ITE=NIT
      IF (IDEFMG.EQ.1) THEN
       MT=MT0
C Michael Koester: If the number is near infinity, don't print it;
C will result in an error
       IF (DEF.LT.1D99) THEN
         RHOLMG=DEF/FD
         RHOLMG=RHOLMG**(1D0/DBLE(NIT))
         IF ((RHOLMG).GE.1D0) BMGEND=.TRUE.

         NITREL = MIN(NIT,IASRLN-1)
         RHOASM=(DEF/RESQUE(1))**(1D0/DBLE(NITREL))

         IF (MT.GE.9) WRITE (6,10000) NIT,DEF,FD,DEF/FD,RHOLMG
         IF (BMSG2) THEN
C           WRITE (CPARAM,'(I15,4D25.16)') NIT,DEF,FD,DEF/FD,RHOLMG
           WRITE (CPARAM,'(I15,3D25.16)') NIT,DEF,DEF/FD,RHOLMG
           CALL OMSG(72,'M011  ')
           WRITE (MPROT,10000) NIT,DEF,FD,DEF/FD,RHOLMG
           WRITE (MPROT,*) 'IER=1 IN M011'
         ENDIF
        ELSE
           WRITE (MPROT,*) 'IER=1 IN M011'
        END IF
      ELSE
       RHOLMG=DBLE(NIT)
      ENDIF
C
      GOTO 99999
C
1000  IER=0
      MT=MT0
      RHOLMG=0D0
      RHOASM=0D0
C Michael Koester: If the number is near infinity, don't print it;
C will result in an error
      IF (DEF.LT.1D99) THEN
        IF (FD.GE.1D-70) RHOLMG=(DEF/FD)**(1D0/DBLE(ITE))
        
        IF (RESQUE(1).GE.1D-70) THEN
          NITREL = MIN(ITE,IASRLN-1)
          RHOASM=(DEF/RESQUE(1))**(1D0/DBLE(NITREL))
        END IF
        
        IF (MT.GE.9) WRITE (6,'(I15,4D25.16)') ITE,DEF,FD,DEF/FD,RHOLMG
        IF (BMSG2) THEN
C         WRITE (CPARAM,'(I15,4D25.16)') ITE,DEF,FD,DEF/FD,RHOLMG
         WRITE (CPARAM,'(I15,3D25.16)') ITE,DEF,DEF/FD,RHOLMG
         IF (DEF.LT.1D99) CALL OMSG(72,'M011  ')
         WRITE (CPARAM,'(D25.16)')  RHOLMG
         IF (RHOLMG.LT.1D99) THEN 
           CALL OMSG(76,'M011  ')
           WRITE (MPROT,10002) ITE,DEF,FD,DEF/FD,RHOLMG
         END IF  
        ENDIF
      ELSE
C DEF=Infinity; RHO=Infinity, set to 1
        RHOLMG = 1D0
        RHOASM = 1D0
        WRITE (MPROT,10001) ITE,DEF
      END IF  
C
99999 MT=MT0 
      IF (BTIME) THEN
       CALL ZTIME(TTMG1)
       TTMG=TTMG+TTMG1-TTMG0
      ENDIF

C Debug: Ausgabe der Anfangsdaten auf feinstem Level
C      print *,'Endvektor:'
C      write (*,'(80D26.18)') (dx(koffx(NLMAX)+i),i=1,80)
C      print *,'Rechte Seite:'
C      write (*,'(80D26.18)') (db(koffb(NLMAX)+i),i=1,80)
C      print *,'Matrix:'
C      write (*,'(80D26.18)') (dwork(L(KLA(NLMAX))+i-1),i=1,80)
C Debug Ende



10000 FORMAT (//' CONVERGENCE OF MG-METHOD FAILED'/
     *        ' NORM OF DEFECT AFTER',I6,' ITERATIONS',D12.3/
     *        ' NORM OF INITIAL DEFECT  ',D12.3/
     *        ' NORM OF DEF/INIT. DEF   ',D12.3/
     *        ' CONVERGENCE RATE        ',D12.3)
10001 FORMAT (' ITERATION',I6,5X,'  RESIDUUM =',D12.3)
10002 FORMAT (' NUMBER OF ITERATIONS    ',I10/
     *        ' NORM OF DEFECT          ',D12.3/
     *        ' NORM OF INITIAL DEFECT  ',D12.3/
     *        ' NORM OF DEF/INIT. DEF   ',D12.3/
     *        ' CONVERGENCE RATE        ',D12.3)
C
C
C
      END
C
