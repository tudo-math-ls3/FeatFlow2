************************************************************************
* FINITE ELEMENT ANALYSIS TOOLBOX  FEAT  (Release 1.0)                 *
*                                                                      *
* Authors: H. Blum, J. Harig, S. Mueller, S. Turek                     *
*          Institute of Applied Mathematics                            *
*          University of Heidelberg                                    *
*          D-6900 HEIDELBERG                                           *
*                                                                      *
************************************************************************
* based on the 2-D routine, modificated by P.Schreiber                 *
* XMOWAn                                                               *
*                                                                      *
* Purpose  Store matrix and pointer vectors corresponding to the       *
*          storage technique  (multigrid version)                      *
*          Successive call of XOWA                                     *
*          Matrix stored in technique  n  (see Reference Manual)       *
*                                                                      *
* Subroutines/functions called  XOWA                                   *
*                                                                      *
* Version from  02/18/91                                               *
*                                                                      *
* INPUT    TYPE                                                        *
* -----    ----                                                        *
* NLEV     I*4    Vector of numbers of the arrays                      *
*                 New arrays are allocated on DWORK for LA(IBLOC)=0    *
* MFILE    I*4    Unit number used for output                          *
* CCFILE   C*(*)  Array of filenames used for output                   *
*                                                                      *
* OUTPUT   TYPE                                                        *
* ------   ----                                                        *
* IER      I*4    Error indicator                                      *
*                 Set by OF0 or OWA                                    *
*                                                                      *
************************************************************************
C
      SUBROUTINE XMOWA3(KLA,KLDIA,KLDIAS,KNDIA,MFILE,CCFILE,IFMT)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER ARR*6
      CHARACTER*(*) CCFILE(*)
      DIMENSION KLA(*),KLDIA(*),KLDIAS(*),KNDIA(*)
C
      PARAMETER (NNARR=299,NNAB=21,NNDER=6,NNLEV=9)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      SAVE /ERRCTL/,/CHAR/,/MGPAR/
C
      SUB='XMOWA3'
      IF (ICHECK.GE.997) CALL OTRC('XMOWA3','02/18/91')
      IER=0
C
      DO 10 ILEV=NLMIN,NLMAX
C
C
      LA   =KLA(ILEV)
      LDIA =KLDIA(ILEV)
      LDIAS=KLDIAS(ILEV)
	NDIA =KNDIA(ILEV)
C
      CALL XOWA(LA   ,'DA    ',MFILE,CCFILE(ILEV),IFMT)
      CALL XOWA(LDIA ,'KDIA  ',MFILE,CCFILE(ILEV),IFMT)
      CALL XOWA(LDIAS,'KDIAS ',MFILE,CCFILE(ILEV),IFMT)
	IF (IFMT.EQ.1) THEN
	 WRITE(MFILE,'(I6)') NDIA
	ELSE
	 WRITE(MFILE) NDIA
	ENDIF
      IF (IER.NE.0) GOTO 99999
      CLOSE(MFILE)
C     
10    CONTINUE   
C     
99999 END
C
C
C
      SUBROUTINE XMOWA7(KLA,KLCOL,KLLD,MFILE,CCFILE,IFMT)
C
      IMPLICIT DOUBLE PRECISION (A,C-H,O-U,W-Z),LOGICAL(B)
      CHARACTER SUB*6,FMT*15,CPARAM*120
      CHARACTER ARR*6
      CHARACTER*(*) CCFILE(*)
      DIMENSION KLA(*),KLCOL(*),KLLD(*)
C
      PARAMETER (NNARR=299,NNAB=21,NNDER=6,NNLEV=9)
      COMMON /ERRCTL/ IER,ICHECK
      COMMON /CHAR/   SUB,FMT(3),CPARAM
      COMMON /MGPAR/  ILEV,NLEV,NLMIN,NLMAX,
     *                ICYCLE,KPRSM(NNLEV),KPOSM(NNLEV)
      SAVE /ERRCTL/,/CHAR/,/MGPAR/
C
      SUB='XMOWA7'
      IF (ICHECK.GE.997) CALL OTRC('XMOWA7','02/18/91')
      IER=0
C
      DO 10 ILEV=NLMIN,NLMAX
C
C
      LA  =KLA(ILEV)
      LCOL=KLCOL(ILEV)
      LLD =KLLD(ILEV)
C
      CALL XOWA(LA  ,'DA    ',MFILE,CCFILE(ILEV),IFMT)
      CALL XOWA(LCOL,'KCOL  ',MFILE,CCFILE(ILEV),IFMT)
      CALL XOWA(LLD ,'KLD   ',MFILE,CCFILE(ILEV),IFMT)
      IF (IER.NE.0) GOTO 99999
      CLOSE(MFILE)
C     
10    CONTINUE   
C     
99999 END
      
      
