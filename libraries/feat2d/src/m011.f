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
* DSTEP    SUBR   DSTEP(DX,DD,DB,DSTEPP)                               *
*                 Returns  DSTEPP := optimal step size for correction  *
* BFILT    BOL    activates the use of filtering techniques            *
*                 by calling DFILT                                     *
* DFILT    SUBR   DFILT(DX,NEQ,IALGP)                                  *
*                 Performs filtering on different positions of the     *
*                 algorithm.                                           *
* KIT0     I*4    auxiliary vectors of length  NLMAX                   *
* KIT      I*4                                                         *
* BNRM     BOL    Check residuum in euclidian norm. If set to .FALSE., *
*                 the norm is measured in the (weighted) l2-norm.      *
* IDEFMG   I*4    Check norm of defect vector for stopping criterion.  *
*                 If =0, no output is written to the console and       *
*                 exactly NIT MG-iterations are performed, regardless  *
*                 of the norm of the defect, i.e. the norm of the      *
*                 defect is not tested for the stopping criterion.     *
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

***********************************************************************
* Description of multigrid filtering by DFILT:
* 
* In:
*  DX     - Solution vector; array [1..NEQ] of double
*  NEQ    - length of solution vector
*  IALGP  - Position in the algorithm
*           0=undefined - enforce the filtering
*           1=on start of the algorithm on finest level; 
*             DX is solution vector
*           2=on start of the MG sweep on the current level;
*             DX is the first solution vector before the
*             sweep starts (0D0-array normally, escept for the
*             finest level, there it's the result of the previous
*             call with IALGP=1) 
*           3=before smoothing; DX is solution vector 
*           4=after smoothing; DX is solution vector 
*           5=before restriction; DX is defect vector
*           6=after restriction; DX is defect vector
*           7=before coarse grid solver;
*             DX is RHS vector on coarsest level
*             (the result of IALGP=6 !)
*           8=before coarse grid solver;
*             DX is start vector on coarsest level (normally filled
*             with 0)
*           9=after coarse grid solver; 
*             DX is calculated solution vector on coarsest level
*          10=before prolongation;
*             DX is update-vector on coarser level
*          11=after prolongation; DX is prolongated update vector
*          12=after coarse grid correction, before post smoothing;
*             DX is solution vector
*          13=after post-smoothing; DX is solution vector
*          14=after one step of MG; DX is the solution vector on the
*             finest level.
*
* Out (to be returned by DFILT):
*  DX     - updated vector.
***********************************************************************

      SUBROUTINE  M011 (DX,DB,DD,KOFFX,KOFFB,KOFFD,KNEQ,NIT,ITE,
     *                  EPS1,EPS2,DEF,DAX,DPROL,DREST,DPRSM,DPOSM,
     *                  DEX,DEXA,DBC,DSTEP,BFILT,DFILT,
     *                  KIT0,KIT,IREL,IDEFMG,RHOLMG,
     *                  RHOASM,BMGEND,BNRM)

      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120

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
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      COMMON /MGTIME/ TTMG,TTS,TTE,TTD,TTP,TTR,IMTIME
      EXTERNAL DAX,DPROL,DREST,DPRSM,DPOSM,DEX,DEXA,DBC,DSTEP
      
      EXTERNAL DFILT

C Debug:
C      PARAMETER (NNARR=299)
C      COMMON          NWORK,IWORK,IWMAX,L(299),DWORK(1)
C      COMMON /MGFLD/  KLA(NNLEV),KLCOLA(NNLEV),
C     *                KLLDA(NNLEV),KLU(NNLEV),
C     *                KLB(NNLEV),KLD(NNLEV),KLAUX(NNLEV)
C      EQUIVALENCE (DWORK(1),VWORK(1),KWORK(1))
C Debug Ende

      SAVE /ERRCTL/,/CHAR/,/OUTPUT/,/MGPAR/,/MGTIME/
      
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

      BMGEND=.FALSE.
      BREL=IREL.EQ.1
      BMSG2=M.GE.2.OR.MT.GE.2
      MT0=MT
      MT=0

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

C *** level counts
C This is due to realise the F/V/W-cycle without recursion

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

C On start of the algorithm perform the first filtering:

      IF (BFILT) CALL DFILT(DX(1+KOFFX(NLMAX)),KNEQ(NLMAX),1)

C Calculate the norm of the defect. The defect vector itself
C will be built in DD:

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
       
C *** FD  is considered as norm of the initial defect

        FD=DEF

C Print out the initial residuum

        IF (BMSG2) THEN
          MT=MT0
          WRITE (CPARAM,'(I15,D25.16)') 0,DEF

C If the number is near infinity, don't print it;
C will result in an error

          IF (DEF.LT.1D99) CALL OMSG(73,'M011  ')

          MT=0
        END IF

C Initialize the queue of the last residuals with DEF

       DO I=1,IASRLN
         RESQUE(I)=DEF
       END DO  

       IF (BTIME) THEN
        CALL ZTIME(TTD1)
        TTD=TTD+TTD1-TTD0
       ENDIF

       IF ((DEF.EQ.0D0).OR.(DEF.LE.EPS2.AND..NOT.BREL)) THEN
C The case "DEF=0" may happen in case of a program error
C (incorrect boundary implementation in the matvec-multiplication).
C This would produce the String "Infinity" in a later call to
C OMSG when printing out the residuum. This is dangerous on some
C machines, as OMSG tries to interpret the string as a number again,
C which obviously fails - and might cause a program crash!
        ITE=0
        GOTO 1000
       ENDIF
      ENDIF

C *** Start multigrid iteration

      DO 100 ITE=1,NIT

C *** initialize level counts for all levels

      DO ILEV=NLMIN,NLMAX
        KIT(ILEV)=KIT0(ILEV)
      END DO

      ILEV=NLMAX
      
C Perform the filtering for the current solution before the MG sweep.

      IF (BFILT) CALL DFILT(DX(1+KOFFX(NLMAX)),KNEQ(NLMAX),2)
      
C On the maximum level we already built out defect vector. If we are
C on a lower level than NLMAX, perform smoothing+restriction down to the
C lowest level NLMIN:
      
110   IF (ILEV.NE.NLMIN) THEN

C First perform filtering of the solution before smoothing.

        IF (BFILT) CALL DFILT(DX(1+KOFFX(ILEV)),KNEQ(ILEV),3)

C Perform the presmoothing with the current solution vector

        IF (BTIME) CALL ZTIME(TTS0)
        IF (KPRSM(ILEV).GT.0) 
     *   CALL DPRSM(DX(1+KOFFX(ILEV)),DB(1+KOFFB(ILEV)),
     *              DD(1+KOFFD(ILEV)),KNEQ(ILEV),KPRSM(ILEV))

        IF (BTIME) THEN
         CALL ZTIME(TTS1)
         TTS=TTS+TTS1-TTS0
        ENDIF

C Perform filtering

        IF (BFILT) CALL DFILT(DX(1+KOFFX(ILEV)),KNEQ(ILEV),4)

C Build the defect vector

        IF (BTIME) CALL ZTIME(TTD0)
        CALL LCP1(DB(1+KOFFB(ILEV)),DD(1+KOFFD(ILEV)),KNEQ(ILEV))
        CALL DAX(DX(1+KOFFX(ILEV)),DD(1+KOFFD(ILEV)),KNEQ(ILEV),
     *           -1D0,1D0)

C Filter the defect vector

        IF (BFILT) CALL DFILT(DD(1+KOFFD(ILEV)),KNEQ(ILEV),5)

        IF (BTIME) THEN
         CALL ZTIME(TTD1)
         TTD=TTD+TTD1-TTD0
        ENDIF

        ILEV=ILEV-1

        IF (BTIME) CALL ZTIME(TTR0)

C Restriction of the defect. The restricted defect is placed
C in DB as the right hand side of the lower level.

        CALL DREST(DD(1+KOFFD(ILEV+1)),DB(1+KOFFB(ILEV)))
       
C eventually skip level ISKPLV

        IF (ILEV.EQ.ISKPLV) THEN
          CALL LCP1(DB(1+KOFFB(ILEV)),DD(1+KOFFD(ILEV)),KNEQ(ILEV))
          ILEV = ILEV -1
          CALL DREST(DD(1+KOFFD(ILEV+1)),DB(1+KOFFB(ILEV)))
        END IF

C Filter the restricted defect vector

        IF (BFILT) CALL DFILT(DB(1+KOFFB(ILEV)),KNEQ(ILEV),6)

        IF (BTIME) THEN
          CALL ZTIME(TTR1)
          TTR=TTR+TTR1-TTR0
        ENDIF

        IF (BTIME) CALL ZTIME(TTD0)

C Choose zero as initial vector on lower level. Implement boundary
C conditions into the just calculated right hand side.

        CALL LCL1(DX(1+KOFFX(ILEV)),KNEQ(ILEV))
        CALL DBC(DB(1+KOFFB(ILEV)),KNEQ(ILEV))
        
C Perform the filtering on the start solution

        IF (BFILT) CALL DFILT(DX(1+KOFFX(ILEV)),KNEQ(ILEV),2)

        IF (BTIME) THEN
          CALL ZTIME(TTD1)
          TTD=TTD+TTD1-TTD0
        ENDIF

C If we are not on the lowest level, repeat the smoothing of the solution/
C restriction of the new defect:

        GOTO 110

      ENDIF

C The previous IF/GOTO sweep ensures that we are on the lowest level now.
C In DD there is the defect on the lowest level, DX is filled with zero
C (plus eventually filtering).

C Do some probably filtering for coarse grid solution and RHS vector

      IF (BFILT) CALL DFILT(DB(1+KOFFB(ILEV)),KNEQ(ILEV),7)
      IF (BFILT) CALL DFILT(DX(1+KOFFX(ILEV)),KNEQ(ILEV),8)

C Solve the system on lowest level:

      IF (BTIME) CALL ZTIME(TTE0)
      CALL DEX(DX(1+KOFFX(NLMIN)),DB(1+KOFFB(NLMIN)),
     *         DD(1+KOFFD(NLMIN)),KNEQ(NLMIN),RHOLMG)

C Filtering the solution on lowest level

      IF (BFILT) CALL DFILT(DX(1+KOFFX(ILEV)),KNEQ(ILEV),9)

      IF (BTIME) THEN
       CALL ZTIME(TTE1)
       TTE=TTE+TTE1-TTE0
      ENDIF

C Prolongate the solution vector, perform the coarse grid correction
C and realise the MG-cycles until we have reached the fine grid again:

130   IF (ILEV.NE.NLMAX) THEN
    
C go to the next higher level
    
        ILEV=ILEV+1

C First perform filtering to the non-prolongated update vector

        IF (BFILT) CALL DFILT(DX(1+KOFFX(ILEV)),KNEQ(ILEV),10)

C Prolongate the update vector; DPROL  returns  DD:=PROL(DX)

        IF (BTIME) CALL ZTIME(TTP0)
        CALL DPROL(DX(1+KOFFX(ILEV-1)),DD(1+KOFFD(ILEV)))

C eventually skip Level ISKPLV

        IF (ILEV.EQ.ISKPLV) THEN
          CALL LCP1 (DD(1+KOFFD(ILEV)),DX(1+KOFFX(ILEV)),KNEQ(ILEV))
          ILEV = ILEV+1
          CALL DPROL(DX(1+KOFFX(ILEV-1)),DD(1+KOFFD(ILEV)))
        END IF

C implement boundary conditions into the prolongated vector
      
        CALL DBC(DD(1+KOFFD(ILEV)),KNEQ(ILEV))

C Perform filtering of the prolongated update vector

        IF (BFILT) CALL DFILT(DD(1+KOFFD(ILEV)),KNEQ(ILEV),11)

C Calculate the step length parameter for the coarse grid correction

        CALL DSTEP(DX(1+KOFFX(ILEV)),DD(1+KOFFD(ILEV)),
     *             DB(1+KOFFB(ILEV)),KNEQ(ILEV),DSTEPP)
     
C Perform the coarse grid correction by adding the coarse grid solution
C (with the calculated step-length parameter) to the current solution
     
        CALL LLC1(DD(1+KOFFD(ILEV)),DX(1+KOFFX(ILEV)),KNEQ(ILEV),
     *            DSTEPP,1D0)

C Perform filtering of the updated solution vector before post-smoothing

        IF (BFILT) CALL DFILT(DX(1+KOFFX(ILEV)),KNEQ(ILEV),12)

        IF (BTIME) THEN
          CALL ZTIME(TTP1)
          TTP=TTP+TTP1-TTP0
        END IF

C *** Post-smoothing

        IF (BTIME) CALL ZTIME(TTS0)
        IF (KPOSM(ILEV).GT.0) 
     *    CALL DPOSM(DX(1+KOFFX(ILEV)),DB(1+KOFFB(ILEV)),
     *               DD(1+KOFFD(ILEV)),KNEQ(ILEV),KPOSM(ILEV))

C Filter the currect solution vector after post-smoothing

        IF (BFILT) CALL DFILT(DX(1+KOFFX(ILEV)),KNEQ(ILEV),13)

        IF (BTIME) THEN
          CALL ZTIME(TTS1)
          TTS=TTS+TTS1-TTS0
        END IF
        
C Update the iteration counter(s) for realising the MG-cycle(s).
C Then either jump to 130 to perform the next prolongation or
C jump to 110 to do perform a next MG sweep on the current
C level.

        KIT(ILEV)=KIT(ILEV)-1
        IF (KIT(ILEV).EQ.0) THEN
          IF (ICYCLE.EQ.0) THEN
            KIT(ILEV)=1
          ELSE
            KIT(ILEV)=KIT0(ILEV)
          END IF
          GOTO 130
        ELSE
          GOTO 110
        END IF

      END IF

CCC      IF (BTIME) CALL ZTIME(TTS0)
CCC      IF (KPRSM(NLMAX).GT.0)  
CCC     * CALL DPRSM(DX(1+KOFFX(NLMAX)),DB(1+KOFFB(NLMAX)),
CCC     *           DD(1+KOFFD(NLMAX)),KNEQ(NLMAX),KPRSM(NLMAX))
CCC
CCC      IF (BTIME) THEN
CCC       CALL ZTIME(TTS1)
CCC       TTS=TTS+TTS1-TTS0
CCC      ENDIF

C We have (hopefully) successfully performed one MG-sweep, starting
C and ending on the finest level. As we are now on the finest level
C again, we can update our defect vector to test the current
C residuum...

C But first perform some possible filtering with the current solution:

      IF (BFILT) CALL DFILT(DX(1+KOFFX(ILEV)),KNEQ(ILEV),14)

      IF (IDEFMG.EQ.1) THEN

        IF (BTIME) CALL ZTIME(TTD0)
       
C Calculate the residuum and its norm; the result can be found in DD:

        CALL LCP1(DB(1+KOFFB(NLMAX)),DD(1+KOFFD(NLMAX)),KNEQ(NLMAX))
        CALL DAX(DX(1+KOFFX(NLMAX)),DD(1+KOFFD(NLMAX)),KNEQ(NLMAX),
     *           -1D0,1D0)

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
C        IF (MT.GT.0) WRITE (6,10001) ITE,DEF

        IF (BMSG2) THEN
          MT=MT0
          WRITE (CPARAM,'(I15,D25.16)') ITE,DEF

C Michael Koester: If the number is near infinity, don't print it;
C will result in an error

          IF (DEF.LT.1D99) CALL OMSG(73,'M011  ')

C        WRITE (MPROT,10001) ITE,DEF
C        WRITE (*,10001) ITE,DEF

          MT=0
        END IF

        IF (BTIME) THEN
          CALL ZTIME(TTD1)
          TTD=TTD+TTD1-TTD0
        ENDIF

C=======================================================================
C ***  Unexpected STOP !!!
C=======================================================================
C Use "not" instead of "GE" because this handles NAN/Infinity-cases
C better!

       IF (.NOT.((DEF.LT.1D99).AND.(DEF.LE.1D3*DEFOLD).AND.
     *    (DEF.LE.5D3*FD))) THEN
        BMGEND=.TRUE.
        GOTO 1000
       ENDIF

      DEFOLD=DEF
       IF (BREL) THEN
        IF (DEF.LE.FD*EPS1.AND.DEF.LE.EPS2.AND.ITE.GE.NIT0) GOTO 1000
       ELSE
        IF (DEF.LE.EPS2) GOTO 1000
       ENDIF
      ENDIF

100   CONTINUE

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

      GOTO 99999

1000  IER=0
      MT=MT0
      RHOLMG=0D0
      RHOASM=0D0
C Michael Koester: If the number is near infinity, don't print it;
C will result in an error
      IF (DEF.LT.1D99) THEN
C FD might be zero in case of incorrect boundary treatment!
C Take care of that!
        IF (FD.EQ.0D0) FD=1D99
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

      END

